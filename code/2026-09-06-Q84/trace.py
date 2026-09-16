"""Session-stable wrapper around tracecheck.emit: one file for the whole PROVE session."""
import sys, os, json
sys.path.insert(0, '/home/clio/projects')
from tracecheck import emit as _E

LOG = '/home/clio/projects/state/trajectory/clio-prove-2026-09-06-Q84.jsonl'

def _attach():
    _E._AGENT = 'clio'
    _E._LOG = LOG
    n = 0
    if os.path.exists(LOG):
        with open(LOG) as f:
            n = sum(1 for line in f if line.strip())
    _E._SEQ = n
    return n

def emit(*a, **k):
    _attach()
    return _E.emit(*a, **k)
