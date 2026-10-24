"""Session-pinned trajectory emitter for PROVE 2026-10-24 Q565.

emit.init()'s default log_dir is RELATIVE and its default session tag is a fresh
timestamp per process.  Both defaults fired today: two files under
projects/state/trajectory/, each restarting ids at evt-0001, so evt-0006's
`inputs: [evt-0001]` resolved to a DIFFERENT event than the one I meant.
This wrapper pins both.  Import it instead of emit.
"""
import sys
sys.path.insert(0, '/home/clio/projects')
from tracecheck.emit import init as _init, emit

LOG_DIR = '/home/clio/state/trajectory'
SESSION = 'Q565-20261024-prove'

def start():
    _init('clio', 'prove', log_dir=LOG_DIR, session=SESSION)

start()
