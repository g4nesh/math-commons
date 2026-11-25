#!/usr/bin/env python3
"""Discovery only. The proof checker never invokes a linear optimizer."""
import sys,json,time,argparse
from pathlib import Path
parser=argparse.ArgumentParser(description='Optional exact phase-I reconstruction of the interior certificate')
parser.add_argument('--output',required=True,type=Path)
args=parser.parse_args()
import coupling as c
F=c.F
keys,A=c.incidence();rr=c.independent_rows(A);u=c.fixture();w0=c.canon(u)
eps=c.FLOOR
b=[sum(x*y for x,y in zip(row,w0))-eps*sum(row)for row in A]
A=[A[i]for i in rr];b=[b[i]for i in rr];sg=[-1 if x<0 else 1 for x in b];n=90;m=len(A)
rows=[[F(si*x)for x in row]+[F(i==j)for j in range(m)]+[si*bb]for i,(row,bb,si)in enumerate(zip(A,b,sg))];basis=list(range(n,n+m));costs=[0]*n+[1]*m
start=time.monotonic()
for loop in range(2000):
 red=[costs[j]-sum(row[j]for row,i in zip(rows,basis)if costs[i])for j in range(n+m)]
 choices=[j for j,v in enumerate(red)if v<0 and j not in basis]
 if not choices:break
 j=min(choices);_,_,i=min((row[-1]/row[j],basis[i],i)for i,row in enumerate(rows)if row[j]>0)
 pivot=[v/rows[i][j]for v in rows[i]]
 for k in range(m):
  if k!=i and rows[k][j]:
   x=rows[k][j];rows[k]=[v-x*y for v,y in zip(rows[k],pivot)]
 rows[i]=pivot;basis[i]=j
 if loop%50==0:print('pivot',loop,'secs',round(time.monotonic()-start,2),flush=True)
else:raise AssertionError('limit')
obj=sum(row[-1]for row,i in zip(rows,basis)if costs[i]);print('final',loop,str(obj),'secs',time.monotonic()-start,flush=True)
if not obj:
 w=[eps]*90
 for row,j in zip(rows,basis):
  if j<90:w[j]+=row[-1]
 assert min(w)>0 and sum(w)==1
 keys,A=c.incidence();assert all(sum(x*y for x,y in zip(row,w))==sum(x*y for x,y in zip(row,w0))for row in A)
 out={'U':[[z.encode()for z in row]for row in u],'colorings':[list(x)for x in c.COLORINGS],'weights':list(map(str,w)),'canonical':list(map(str,w0)),'minimumWeight':str(min(w)),'correction':[str(x-y)for x,y in zip(w,w0)]}
 args.output.write_text(json.dumps(out,indent=2)+'\n')
