"""Anonymous Orion range bridge with durable pre-I/O charges and pinned files."""
import runtime as rt
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import re
import threading
import urllib.parse

import requests

REPOSITORY = 'slaf-project/X-Atlas-Orion'


class PublicStore:
    def __init__(self, budget, session=None):
        self.budget = budget
        self.http = session or requests.Session()
        self.http.trust_env = False
        self.links = {}
        self.lock = threading.Lock()
        self.errors = 0
        self.requests = 0
        revision = budget.bindings.get('repository_revision')
        if revision is None:
            # Bound even failed API requests before opening the connection.
            budget.charge(1_000_000, 'repository revision metadata upper bound')
            with self.http.get('https://huggingface.co/api/datasets/' + REPOSITORY,
                               headers={'Accept-Encoding': 'identity'},
                               stream=True, allow_redirects=False, timeout=(10, 60)) as response:
                response.raise_for_status()
                # raw.read bypasses requests' normal content decoding. Decode
                # even if an upstream returns gzip despite requesting identity,
                # and bound the decoded JSON instead of using response.json().
                if 'Content-Length' in response.headers:
                    assert int(response.headers['Content-Length']) < 1_000_000
                body = response.raw.read(1_000_000, decode_content=True)
                assert len(body) < 1_000_000, 'Repository metadata exceeds reserve'
                revision = rt.json.loads(body)['sha']
            assert re.fullmatch(r'[0-9a-f]{40}', revision)
            budget.bind('repository_revision', revision)
        assert re.fullmatch(r'[0-9a-f]{40}', revision)
        self.revision = revision
        parent = self

        class Handler(BaseHTTPRequestHandler):
            protocol_version = 'HTTP/1.1'

            def log_message(self, *args):
                pass

            def do_HEAD(self):
                self.serve(True)

            def do_GET(self):
                self.serve(False)

            def serve(self, head):
                sent_headers = False
                try:
                    path = urllib.parse.unquote(urllib.parse.urlsplit(self.path).path)
                    assert path.startswith('/public/')
                    key = path[len('/public/'):]
                    assert '..' not in key and key.startswith('data/')
                    url, source = parent.resolve(key)
                    size = source['size']
                    if head:
                        self.send_response(200)
                        self.send_header('Content-Length', str(size))
                        self.send_header('ETag', source['etag'])
                        self.end_headers()
                        return
                    start, end, partial = parent.bounds(self.headers.get('Range'), size)
                    amount = end - start + 1
                    rt.guard()
                    parent.budget.charge(amount, 'Orion GET attempt ' + key)
                    headers = {'Accept-Encoding': 'identity', 'If-Match': source['etag']}
                    if partial:
                        headers['Range'] = f'bytes={start}-{end}'
                    with parent.http.get(url, headers=headers, stream=True,
                                         allow_redirects=False, timeout=(10, 90)) as response:
                        assert response.status_code == (206 if partial else 200)
                        assert int(response.headers['Content-Length']) == amount
                        assert response.headers['ETag'] == source['etag']
                        assert response.headers.get('Content-Encoding', 'identity') == 'identity'
                        if partial:
                            assert response.headers['Content-Range'] == f'bytes {start}-{end}/{size}'
                        self.send_response(response.status_code)
                        self.send_header('Content-Length', str(amount))
                        self.send_header('ETag', source['etag'])
                        if partial:
                            self.send_header('Content-Range', response.headers['Content-Range'])
                        self.end_headers()
                        sent_headers = True
                        received = 0
                        while received < amount:
                            block = response.raw.read(min(1024 * 1024, amount - received))
                            assert block, 'Truncated Orion payload'
                            received += len(block)
                            self.wfile.write(block)
                        with parent.lock:
                            parent.requests += 1
                except Exception:
                    with parent.lock:
                        parent.errors += 1
                    self.close_connection = True
                    if not sent_headers:
                        try:
                            self.send_error(503, 'Public range refused')
                        except (BrokenPipeError, ConnectionResetError):
                            pass

        self.server = ThreadingHTTPServer(('127.0.0.1', 0), Handler)
        self.thread = threading.Thread(target=self.server.serve_forever, daemon=True)
        self.thread.start()

    @staticmethod
    def bounds(value, size):
        assert size > 0
        if value is None:
            return 0, size - 1, False
        match = re.fullmatch(r'bytes=(\d*)-(\d*)', value)
        assert match and any(match.groups()), 'Unsupported byte range'
        first, last = match.groups()
        if first:
            start = int(first)
            end = min(int(last), size - 1) if last else size - 1
        else:
            assert int(last) > 0
            start, end = max(0, size - int(last)), size - 1
        assert 0 <= start <= end < size
        return start, end, True

    def resolve(self, key):
        with self.lock:
            if key not in self.links:
                self.budget.charge(64_000, 'Orion HEAD metadata upper bound ' + key)
                url = 'https://huggingface.co/datasets/' + REPOSITORY + '/resolve/' + self.revision + '/' + key
                with self.http.head(url, allow_redirects=True, timeout=(10, 60)) as response:
                    response.raise_for_status()
                    assert response.status_code == 200
                    source = dict(size=int(response.headers['Content-Length']), etag=response.headers['ETag'])
                    assert source['size'] > 0 and source['etag']
                    self.budget.bind('file:' + key, source)
                    self.links[key] = response.url, source
            return self.links[key]

    def options(self):
        return dict(endpoint=f'http://127.0.0.1:{self.server.server_port}',
                    allow_http='true', skip_signature='true', region='us-east-1',
                    virtual_hosted_style_request='false')

    def close(self):
        self.server.shutdown()
        self.server.server_close()
        self.thread.join()
        self.http.close()
