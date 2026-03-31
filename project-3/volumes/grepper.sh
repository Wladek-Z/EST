#!/bin/bash
grep "Final energy" iceIh-V/*.out > iceIh-energies.txt
grep "unit-cell volume" iceIh-V/*.out > iceIh-volumes.txt
grep "Final energy" iceII-V/*.out > iceII-energies.txt
grep "unit-cell volume" iceII-V/*.out > iceII-volumes.txt
grep "Final energy" iceVIII-V/*.out > iceVIII-energies.txt
grep "unit-cell volume" iceVIII-V/*.out > iceVIII-volumes.txt
