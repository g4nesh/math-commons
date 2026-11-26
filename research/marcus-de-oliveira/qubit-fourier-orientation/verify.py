"""Replay a sufficient all-axis polytope certificate without numerical libraries."""
import json
from fractions import Fraction as F
from pathlib import Path
from model import FACETS,ORBIT_REPS,ROWS,orbit_column,require,target,vertices


def verify(data):
    require(data.get('schema')=='qubit-fourier-orientation-v1','certificate schema')
    expected=set(vertices())
    supplied={}
    columns=tuple(orbit_column(c) for c in ORBIT_REPS)
    require(len(ORBIT_REPS)==30 and len(ROWS)==358,'model dimensions')
    for item in data['vertices']:
        require(all(isinstance(x,str) for x in item['point']+item['weights']),'certificate scalars must be exact strings')
        v=tuple(F(x) for x in item['point']);w=tuple(F(x) for x in item['weights'])
        require(len(v)==4 and len(w)==30,'certificate dimensions')
        require(v not in supplied,'duplicate vertex')
        require(v in expected,'unexpected vertex')
        require(all(x>=0 for x in w) and sum(w)==1,'probability weights')
        b=target(v)
        require(all(sum(col[j]*mass for col,mass in zip(columns,w))==b[j] for j in range(len(b))),
                'histogram identity')
        supplied[v]=w
    require(set(supplied)==expected,'missing vertex')
    return supplied


def null_direction(active):
    """A nonzero rational vector orthogonal to all active facet normals."""
    a=[[F(x) for x in row] for row in active];pivots=[];row=0
    for j in range(4):
        i=next((i for i in range(row,len(a)) if a[i][j]),None)
        if i is None:continue
        a[row],a[i]=a[i],a[row];z=a[row][j];a[row]=[x/z for x in a[row]]
        for i in range(len(a)):
            if i!=row:
                z=a[i][j];a[i]=[x-z*y for x,y in zip(a[i],a[row])]
        pivots.append(j);row+=1
    free=next((j for j in range(4) if j not in pivots),None)
    if free is None:return None
    result=[F(0)]*4;result[free]=F(1)
    for i,j in enumerate(pivots):result[j]=-a[i][free]
    return tuple(result)


def construct(point,vertex_laws):
    """Exact rational query by face-preserving line splits (at most4 levels).

    The same arithmetic algorithm is valid over real inputs; this executable
    deliberately accepts rational strings/integers/Fractions, not float guesses.
    """
    require(len(point)==4 and all(isinstance(x,(str,int,F)) and not isinstance(x,bool) for x in point),
            'use exact rational coordinates')
    p=tuple(F(x) for x in point)
    def recurse(p,depth):
        slack=[b-sum(x*y for x,y in zip(a,p)) for a,b in FACETS]
        require(all(x>=0 for x in slack),'point outside sufficient polytope')
        if p in vertex_laws:return vertex_laws[p]
        require(depth<4,'face dimension did not decrease')
        direction=null_direction([a for (a,_),s in zip(FACETS,slack) if not s])
        require(direction is not None,'missing full-rank vertex')
        speeds=[sum(x*y for x,y in zip(a,direction)) for a,_ in FACETS]
        positive=min(s/v for s,v in zip(slack,speeds) if v>0)
        negative=min(-s/v for s,v in zip(slack,speeds) if v<0)
        require(positive>0 and negative>0,'non-interior line in current face')
        plus=tuple(x+positive*y for x,y in zip(p,direction))
        minus=tuple(x-negative*y for x,y in zip(p,direction))
        a=recurse(plus,depth+1);b=recurse(minus,depth+1)
        return tuple((negative*x+positive*y)/(positive+negative) for x,y in zip(a,b))
    return recurse(p,0)


if __name__=='__main__':
    data=json.loads(Path(__file__).with_name('certificate.json').read_text())
    laws=verify(data)
    for p in [(0,0,0,0),(0,0,0,F(7,4)),(1,1,1,0),(F(1,5),F(2,7),F(-1,8),F(1,2))]:
        w=construct(p,laws);b=target(p)
        require(all(x>=0 for x in w) and sum(w)==1,'interpolated law')
        columns=[orbit_column(c) for c in ORBIT_REPS]
        require(all(sum(col[j]*mass for col,mass in zip(columns,w))==b[j] for j in range(len(b))),
                'interpolated histogram')
    print(f'PASS: {len(laws)} complete vertices, 30 orbits, 358 exact coefficients per vertex; 4 exact interpolation queries.')
