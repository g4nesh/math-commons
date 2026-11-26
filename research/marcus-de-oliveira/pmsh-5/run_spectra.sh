#!/bin/sh
# Special n = 5 spectra (logs/n5_special.txt). Usage: PY=python3 ./run_spectra.sh
for spec in 0,0.25,0.5,0.75,1 0,0.01,0.02,0.03,1 0,0.49,0.5,0.51,1 0,0.001,0.01,0.1,1 0,0.01,0.5,0.99,1 0,0.3,0.3001,0.3002,1 0,0.0001,0.5,0.9999,1 0,0.5,0.5,0.75,1 ; do
  "$PY" facet_search.py 5 7 1500 6 $spec
done
