"""Compare cached reductions with unmodified official metric/anchor outputs."""
import runtime as rt
import runpy,numpy as np
import stat_cache,bounded_score,official_scoring
stat_cache.install_unchecked();genes=np.array(['a','b','c']);x=np.array([[1,4,2],[3,2,7],[4,1,2]],np.int32);cfg=official_scoring.config();left=bounded_score.make(x,'first',genes);right=bounded_score.make(x,'second',genes);stat_cache.cached(left,cfg);actual=stat_cache.cached(right,cfg);expected=stat_cache.ORIGINAL(right,cfg)
for a,b in zip(actual,expected):np.testing.assert_array_equal(a,b)
runpy.run_path('verify_bounded_score.py',run_name='__main__');assert stat_cache.hits>0 and stat_cache.misses>0;rt.save('artifacts/stat_cache_fixture.json',dict(status='PASS bitwise metrics, anchors and relabeled single-group moments',implementation_sha256=rt.sha('stat_cache.py'),official_files=official_scoring.provenance()['official_files'],hits=stat_cache.hits,misses=stat_cache.misses,scope='Existing full official fixture plus explicit relabeling invariance; stored numeric outputs are lossless'));print('PASS exact official stat cache',stat_cache.hits,'hits',stat_cache.misses,'misses')
