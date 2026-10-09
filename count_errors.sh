#!/bin/bash
# awk counter that ALWAYS prints an integer (grep -c exits 1 on zero matches and
# command substitution swallows it -- blank and 0 are the same pixel, different facts)
awk 'BEGIN{n=0} /^!/{n++} END{printf "%d\n", n}' "$1"
