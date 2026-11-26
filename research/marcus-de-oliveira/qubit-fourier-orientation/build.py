"""Optional numerical discovery; emitted candidates require verify.py exact replay."""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from fractions import Fraction as F
from model import ORBIT_REPS,orbit_column,target,vertices,require

columns=tuple(orbit_column(c) for c in ORBIT_REPS)
A=np.array(columns,float).T
out=[]
for v in vertices():
    b=target(v)
    result=linprog(np.zeros(len(columns)),A_eq=A,b_eq=np.array(b,float),bounds=(0,None),method='highs')
    require(result.success,'LP candidate missing')
    weights=tuple(F(float(x)).limit_denominator(1000000) for x in result.x)
    require(all(x>=0 for x in weights),'negative candidate')
    require(all(sum(col[j]*w for col,w in zip(columns,weights))==b[j] for j in range(len(b))),'inexact candidate')
    out.append({'point':[str(x) for x in v],'weights':[str(x) for x in weights]})
p=Path(__file__).with_name('certificate.json')
p.write_text(json.dumps({'schema':'qubit-fourier-orientation-v1','vertices':out},indent=2)+'\n')
print('Exact candidate lifts saved:',len(out))
