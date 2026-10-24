"""Stable trace handle for this session: absolute log_dir + pinned session tag,
so every process appends to ONE file and evt-ids continue instead of colliding.
(emit.py's own docstring warns about both faults; both have happened here before.)"""
import sys
sys.path.insert(0, '/home/clio/projects')
from tracecheck.emit import init, emit
LOG = init("clio", "prove",
           log_dir="/home/clio/projects/state/trajectory",
           session="Q565-20261024c2")
