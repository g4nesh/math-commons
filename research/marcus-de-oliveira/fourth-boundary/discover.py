#!/usr/bin/env python3
"""Optional exact phase-I discovery; the main certificate replay needs no LP."""
from fractions import Fraction as F
from itertools import permutations
from pathlib import Path
import json
import sys
import boundary
m = boundary.m
P = boundary.P
def system(t):return m.coefficient_system(m.s.Q(t,-t),P)
def phase(t):
 A,b=system(t); labels=list(range(len(A)));signs=[-1 if v<0 else 1 for v in b]; n=24;nr=len(A)
 rows=[[sg*x for x in r]+[F(i==j)for j in range(nr)]+[sg*y] for i,(r,y,sg)in enumerate(zip(A,b,signs))]
 basis=list(range(n,n+nr));costs=[0]*n+[1]*nr
 for loop in range(2000):
  red=[costs[j]-sum(row[j] for row,i in zip(rows,basis)if costs[i])for j in range(n+nr)]
  choices=[j for j,v in enumerate(red)if v<0 and j not in basis]
  if not choices:break
  j=min(choices);choices=[(row[-1]/row[j],basis[i],i)for i,row in enumerate(rows)if row[j]>0];_,_,i=min(choices)
  pivot=[v/rows[i][j]for v in rows[i]]
  for k in range(nr):
   if k!=i and rows[k][j]:
    ratio=rows[k][j];rows[k]=[v-ratio*w for v,w in zip(rows[k],pivot)]
  rows[i]=pivot;basis[i]=j
 else:raise RuntimeError('limit')
 obj=sum(row[-1]for row,i in zip(rows,basis)if costs[i])
 w=[F(0)]*24
 for row,i in zip(rows,basis):
  if i<24:w[i]=row[-1]
 dual=[sg*sum(row[n+j]for row,i in zip(rows,basis)if costs[i])for j,sg in enumerate(signs)]
 assert all(sum(y*r[j]for y,r in zip(dual,A))<=0 for j in range(24))
 assert sum(y*v for y,v in zip(dual,b))==obj
 if not obj:assert all(sum(x*y for x,y in zip(r,w))==v for r,v in zip(A,b))
 return {'s':str(t),'objective':str(obj),'support':[list(P[j])for j,wj in enumerate(w)if wj], 'weights':[str(wj)for wj in w if wj],'dual':list(map(str,dual)),'basis':basis,'loops':loop}
if __name__=='__main__':
 vals=list(map(F,sys.argv[1:]))or[F(i,20)for i in range(21)]
 for t in vals:
  r=phase(t);print(json.dumps(r),flush=True)
