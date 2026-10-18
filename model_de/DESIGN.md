# DE experiments

strategy.py implements source-sign ranking, fixed/adaptive support masks, selective
control templates and integer redistribution. exchange.py performs exchanges within
equal-depth groups. exchange_verify.py checks integer row/column margins independently.

harness.py evaluates destination-family holdouts with 400 predicted cells per target
and the official six metrics. It needs user-prepared public tasks and source banks.
Do not use held-out treated labels for fitting, support selection or predictions.
Local DE call counts measure predictions rather than their truth accuracy.

The optional --fast rank kernel is excluded pending an upstream license. Use the
default official backend. verify.py retains synthetic count, metric and parity tests.
