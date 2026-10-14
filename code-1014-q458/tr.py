"""Stable-tag trajectory emitter for this session. Absolute log_dir + pinned session,
per the two faults recorded in tracecheck/emit.py's own docstring."""
import sys
sys.path.insert(0, '/home/clio/projects')
from tracecheck.emit import init, emit as _emit
init("clio", "prove", log_dir="/home/clio/projects/state/trajectory",
     session="Q458-20261014c2")
def emit(*a, **k): return _emit(*a, **k)
