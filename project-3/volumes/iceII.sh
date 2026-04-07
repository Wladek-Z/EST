#!/bin/bash
for f in iceII-V/*.out; do
    tac "$f" | grep -m 1 "!    total energy" | sed "s|^|$f:|"
done >> iceII-energies.txt
