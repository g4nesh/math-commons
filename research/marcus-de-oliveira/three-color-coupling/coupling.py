"""Exact finite coloring statistics and canonical signed weights."""
from fractions import Fraction as F
from itertools import combinations,permutations,product
from pathlib import Path
import importlib.util,sys,json
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('color_exact',HERE.parent/'certificates'/'verify.py');e=importlib.util.module_from_spec(spec);sys.modules[spec.name]=e;spec.loader.exec_module(e)
Q=e.QComplex
Z=Q(0);O=Q(1)
PATTERN=(2,2,2);BLOCKS=((0,1),(2,3),(4,5));N=6
COLORINGS=tuple(x for x in product(range(3),repeat=N)if tuple(x.count(j)for j in range(3))==PATTERN)
SUBSETS=tuple(s for k in range(N+1)for s in combinations(range(N),k))

def det(a):return e.determinant(a)
def minor(u,rows,cols):return det([[u[i][j]for j in cols]for i in rows])
def canon(u):
 out=[];du=det(u)
 for coloring in COLORINGS:
  groups=[[i for i,c in enumerate(coloring)if c==j]for j in range(3)]
  order=sum(groups,[]);sign=(-1)**sum(order[i]>order[j]for i in range(N)for j in range(i+1,N))
  z=Q(sign)*du.conjugate()
  for rows,cols in zip(groups,BLOCKS):z=z*minor(u,rows,cols)
  out.append(z.real)
 assert sum(out)==1
 return out

def incidence():
 keys=[];rows=[]
 for subset in SUBSETS:
  counts=sorted(set(tuple(sum(c[i]==j for i in subset)for j in range(3))for c in COLORINGS))
  for count in counts:
   keys.append((subset,count));rows.append([int(tuple(sum(c[i]==j for i in subset)for j in range(3))==count)for c in COLORINGS])
 return keys,rows

def independent_rows(rows):
 pivots={};selected=[]
 for i,row in enumerate(rows):
  v=list(map(F,row))
  for col,p in pivots.items():
   if v[col]:
    a=v[col];v=[x-a*y for x,y in zip(v,p)]
  col=next((j for j,x in enumerate(v)if x),None)
  if col is not None:
   factor=v[col];pivots[col]=[x/factor for x in v];selected.append(i)
 return selected

def rotation(n,i,j,t,phase=O):
 c=(1-t*t)/(1+t*t);ss=2*t/(1+t*t)
 g=e.identity(n);g[i][i]=g[j][j]=Q(c);g[i][j]=Q(ss)*phase;g[j][i]=-Q(ss)*phase.conjugate();return g

def fixture():
 # Embedded archived four-cycle, followed by explicit small complex mixing.
 nums=[[-31,492,-240,300],[-300,225,500,0],[240,-180,225,500],[-492,-256,-180,225]]
 u=e.identity(6)
 for i in range(4):
  for j in range(4):u[i][j]=Q(F(nums[i][j],625))
 for k,(i,j)in enumerate([(0,4),(1,5),(2,4),(3,5),(0,5),(1,4),(2,5),(3,4),(4,5)]):
  u=e.multiply(rotation(6,i,j,F(1,100),Q(0,1)if k%2 else O),u)
 return u


COVER = [((0,),(0,)),((1,),(3,)),((2,),(1,)),((3,),(2,))] + [
    (r,t) for r,ts in [((0,1),[(0,2),(0,3),(1,3)]),
                       ((0,2),[(0,1),(0,3),(1,2)]),
                       ((0,3),[(0,1),(0,2),(2,3)])] for t in ts]
COVER += [((i,),(j,)) for i in (4,5) for j in range(4)]
FLOOR=F(1,10**24)

def abs2(z):return z.real*z.real+z.imag*z.imag

def cover_sum(u):return sum(abs2(minor(u,r,t))for r,t in COVER)

def check_unitary(u):
    if len(u)!=6 or any(len(row)!=6 or any(not isinstance(z,Q)for z in row)for row in u):
        raise ValueError('Expected an exact six by six Gaussian-rational matrix')
    if e.multiply(u,e.adjoint(u))!=e.identity(6) or e.multiply(e.adjoint(u),u)!=e.identity(6):
        raise ValueError('Matrix is not exactly unitary')

def load_data():return json.loads((HERE/'interior.json').read_text())

def corrected_weights(u,data=None):
    data=load_data()if data is None else data
    verify(data)
    check_unitary(u);base=fixture();gap=1-cover_sum(base);radius=min(FLOOR/60,gap/192)
    if any(abs2(x-y)>radius*radius for r,t in zip(u,base)for x,y in zip(r,t)):
        raise ValueError('Matrix lies outside the certified entrywise ball')
    correction=list(map(F,data['correction']));values=[x+y for x,y in zip(canon(u),correction)]
    if len(values)!=90 or min(values)<FLOOR/2 or sum(values)!=1:
        raise ValueError('Corrected weights failed their certified bounds')
    return values

def verify(data=None):
    data=load_data()if data is None else data
    u=fixture();check_unitary(u)
    if data['U']!=[[z.encode()for z in row]for row in u]:raise ValueError('Archived fixture changed')
    if data['colorings']!=[list(x)for x in COLORINGS]:raise ValueError('Wrong fixed-count coloring list')
    for key in ('weights','canonical','correction'):
        if len(data[key])!=90 or any(type(x)is not str for x in data[key]):raise ValueError('Expected90 exact rational strings')
    w=list(map(F,data['weights']));canonical=canon(u);correction=list(map(F,data['correction']))
    if min(w)!=FLOOR or sum(w)!=1:raise ValueError('Base convex weights have wrong floor or sum')
    if canonical!=list(map(F,data['canonical'])):raise ValueError('Canonical determinant weights changed')
    if correction!=[x-y for x,y in zip(w,canonical)]:raise ValueError('Correction coefficients changed')
    keys,rows=incidence()
    if any(sum(x*y for x,y in zip(row,correction))!=0 for row in rows):raise ValueError('Correction is not in the grouped-statistic kernel')
    # Independent route: all 924 individual squared minors, grouped by colors.
    grouped={key:F(0)for key in keys};minor_count=0
    for subset in SUBSETS:
        for columns in combinations(range(6),len(subset)):
            count=tuple(sum(j in block for j in columns)for block in BLOCKS)
            grouped[(subset,count)]+=abs2(minor(u,subset,columns));minor_count+=1
    for key,row in zip(keys,rows):
        if sum(x*y for x,y in zip(row,w))!=grouped[key]:raise ValueError('Direct grouped squared-minor identity failed')
    scores=[sum(tuple(sorted(p[j]for j in t))==r for r,t in COVER)for p in permutations(range(6))]
    if min(scores)!=1:raise ValueError('The cover fails on a permutation')
    gap=1-cover_sum(u)
    if gap<=0:raise ValueError('No strict common-minor obstruction')
    radius=min(FLOOR/60,gap/192)
    if FLOOR-30*radius<FLOOR/2 or gap-96*radius<gap/2:raise ValueError('Joint neighborhood bounds failed')
    if data['minimumWeight']!=str(FLOOR):raise ValueError('Recorded floor mismatch')
    return {'dimension':6,'pattern':[2,2,2],'colorings':90,'groupedStatisticRows':len(keys),
            'individualSquaredMinors':minor_count,'strictBaseMassFloor':str(FLOOR),
            'negativeCanonicalMasses':sum(x<0 for x in canonical),'zeroMatrixEntries':sum(z==Z for row in u for z in row),
            'permutationsChecked':len(scores),'coverTerms':len(COVER),'minimumPermutationCover':min(scores),
            'baseObstructionGap':str(gap),'closedEntrywiseUnitaryBallRadius':str(radius),
            'massLowerBoundThroughoutBall':str(FLOOR/2),'obstructionGapLowerBoundThroughoutBall':str(gap/2),
            'massLipschitzConstant':30,'coverLipschitzConstant':96,
            'scope':'Explicit nonempty open subclass of pattern(2,2,2), outside the common-minor polytope throughout. Not general Conjecture F.'}


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report',type=Path)
    args=parser.parse_args();report=verify()
    if args.report:args.report.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
