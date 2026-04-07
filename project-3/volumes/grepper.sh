#!/bin/bash
grep "unit-cell volume" iceIh-V/*.out > iceIh-volumes.txt
grep "unit-cell volume" iceII-V/*.out > iceII-volumes.txt
grep "unit-cell volume" iceVIII-V/*.out > iceVIII-volumes.txt
