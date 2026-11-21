#!/bin/bash
cd /home/clio/projects/proofs/code-1121c2-Q597
{ time python3 shadowlaw.py 5 4 12 > shadowlaw5.out ; } 2> shadowlaw5.time
{ time KEEP=3 python3 shadowlaw.py 6 4 24 > shadowlaw6.out ; } 2> shadowlaw6.time
