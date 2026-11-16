"""Minimal runner (no pytest needed): runs every test_* function in tests/test_brdfkit.py."""
import sys
import time
import traceback

import test_brdfkit as T

names = [n for n in dir(T) if n.startswith('test_')]
failed = 0
for n in names:
    t0 = time.time()
    try:
        r = getattr(T, n)()
        print(f'  ok    {n:32s} {time.time()-t0:5.1f}s  {r or ""}')
    except Exception:
        failed += 1
        print(f'  FAIL  {n}'); traceback.print_exc()
print(f'{len(names) - failed} passed, {failed} failed')
sys.exit(1 if failed else 0)
