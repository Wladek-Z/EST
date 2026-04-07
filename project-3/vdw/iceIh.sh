#!/bin/bash
for f in iceIh-V/*.out; do
    tac "$f" | grep -m 1 "!    total energy" | sed "s|^|$f:|"
done >> iceIh-energies.txt
