#!/bin/bash
cd /home/clio/projects/proofs/code-1121c2-Q597
{ time KEEP=3 python3 handshake.py 6 4 24 > handshake6.out ; } 2> handshake6.time
