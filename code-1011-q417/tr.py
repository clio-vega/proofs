"""Trajectory emitter, with an ABSOLUTE log_dir and a STABLE session tag.
The bare emit.init() default is relative (see emit.py's capitalised warning) and
scattered one trajectory across three directories on 2026-10-09."""
import sys
sys.path.insert(0, '/home/clio/projects')
from tracecheck.emit import init as _init, emit

def init():
    _init("clio", "prove",
          log_dir="/home/clio/projects/state/trajectory",
          session="Q417-20261011c1")
