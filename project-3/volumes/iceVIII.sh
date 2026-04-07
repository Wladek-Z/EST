#!/bin/bash
for f in iceVIII-V/*.out; do
    tac "$f" | grep -m 1 "!    total energy" | sed "s|^|$f:|"
done >> iceVIII-energies.txt
