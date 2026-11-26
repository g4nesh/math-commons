"""Orbit A on U(5): exact gap decomposition and equality mechanisms.

Python 3.9+ standard library; all mathematical arithmetic is rational.
"""
from fractions import Fraction as F
from itertools import combinations, permutations
from pathlib import Path
import importlib.util
import json
import sys

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('physical_a_exact',HERE.parents[2]/'certificates'/'verify.py')
exact=importlib.util.module_from_spec(spec)
sys.modules[spec.name]=exact
spec.loader.exec_module(exact)
Q=exact.QComplex
PERMS=tuple(permutations(range(5)))
LOOKUP={p:i for i,p in enumerate(PERMS)}
PAIRS=tuple((i,j)for k in (1,2)for i in combinations(range(5),k)for j in combinations(range(5),k))


def require(ok,message):
    if not ok:raise ValueError(message)


def squared(z):return z.real*z.real+z.imag*z.imag


def minor(u,rows,cols):return exact.determinant([[u[i][j]for j in cols]for i in rows])


def gap_data(u):
    require(len(u)==5 and all(len(row)==5 and all(isinstance(x,Q)for x in row)for row in u),'exact 5 by 5 Gaussian matrix')
    require(exact.multiply(u,exact.adjoint(u))==exact.identity(5),'exact unitarity')
    a=u[0][0];p=minor(u,(0,3),(3,4));q=minor(u,(0,4),(3,4))
    y=[a,p,q]
    c=[minor(u,(3,4),(3,4)),u[4][0],-u[3][0]]
    det=minor(u,(0,3,4),(0,3,4))
    require(det==sum((v*w for v,w in zip(y,c)),Q(0)),'Laplace identity')
    require(squared(det)==squared(minor(u,(1,2),(1,2))),'complementary unitary minors')
    deficit_minors=[minor(u,(3,4),cols)for cols in combinations((1,2,3,4),2)if cols!=(3,4)]
    tau=sum(squared(x)for x in deficit_minors)
    require(tau==1-sum(squared(x)for x in c),'two-row Cauchy-Binet deficit')
    crosses=[y[i]*c[j].conjugate()-y[j]*c[i].conjugate()for i,j in combinations(range(3),2)]
    norm_y=sum(squared(x)for x in y)
    gap=norm_y-squared(det)
    require(gap==tau*norm_y+sum(squared(x)for x in crosses),'exact sum-of-squares identity')
    require(gap>=0,'physical orbit-A inequality')
    zero_rhs=all(x==Q(0)for x in y)
    block_zero=all(u[i][j]==Q(0)for i in (3,4)for j in (1,2))
    column_supported=all(u[i][0]==Q(0)for i in (0,1,2))
    require((tau==0)==(block_zero or column_supported),'deficit mechanism classification')
    require((gap==0)==(zero_rhs or block_zero or column_supported),'complete physical equality classification')
    # Connected row-overlap graph witnesses irreducibility under permutation block form.
    seen={0}
    while True:
        added={j for i in seen for j in range(5)if any(u[i][k]!=Q(0)and u[j][k]!=Q(0)for k in range(5))}
        if added<=seen:break
        seen|=added
    return {'gap':str(gap),'tau':str(tau),'rhs':str(norm_y),
            'zero_rhs':zero_rhs,'block_zero':block_zero,'column_supported':column_supported,
            'connected_row_overlap':len(seen)==5,
            'cross_squared':list(map(str,(squared(x)for x in crosses)))}


def householder(vector):
    norm=sum(squared(x)for x in vector)
    require(norm>0,'nonzero Householder vector')
    return [[Q(int(i==j))-2*vector[i]*vector[j].conjugate()/Q(norm)for j in range(len(vector))]for i in range(len(vector))]


def embed(block,indices):
    result=exact.identity(5)
    for i,row in enumerate(block):
        for j,x in enumerate(row):result[indices[i]][indices[j]]=x
    return result


def fixtures():
    block=exact.multiply(embed(householder([Q(1),Q(1,1),Q(2,-1)]),(0,1,2)),
                         embed(householder([Q(1,1),Q(2,-1),Q(1)]),(0,3,4)))
    v4=householder([Q(1,1),Q(2,-1),Q(-1,2),Q(3,1)])
    supported=[[Q(0)for _ in range(5)]for _ in range(5)]
    supported[3][0]=Q(1)
    for i,row in zip((0,1,2,4),v4):
        for j,x in zip((1,2,3,4),row):supported[i][j]=x
    gate=[[Q(F(3,5)),Q(0,F(4,5))],[Q(0,F(4,5)),Q(F(3,5))]]
    supported=exact.multiply(embed(gate,(3,4)),supported)
    swap=exact.identity(5);swap[0],swap[1]=swap[1],swap[0]
    zero=exact.multiply(exact.multiply(embed(v4,(1,2,3,4)),swap),embed(gate,(1,2)))
    dense=exact.identity(5)
    for v in ([Q(1,1),Q(2,-1),Q(-1,2),Q(3),Q(0,1)],
              [Q(2),Q(-1,1),Q(1,-2),Q(0,2),Q(3,1)]):
        dense=exact.multiply(householder(v),dense)
    return {'block_zero_contact':block,'column_supported_contact':supported,'zero_rhs_contact':zero,'strict_dense':dense}


def indicator(p,k,j):return int(tuple(sorted(p[i]for i in k))==j)


def gap_permutation(p):
    return (indicator(p,(0,),(0,))+indicator(p,(0,3),(3,4))+
            indicator(p,(0,4),(3,4))-indicator(p,(1,2),(1,2)))


def rank(rows):
    pivots={}
    for source in rows:
        row=list(map(F,source))
        for j,pivot in pivots.items():
            if row[j]:
                value=row[j];row=[x-value*y for x,y in zip(row,pivot)]
        j=next((j for j,x in enumerate(row)if x),None)
        if j is not None:
            value=row[j];pivots[j]=[x/value for x in row]
    return len(pivots)


def facet_data():
    values=[gap_permutation(p)for p in PERMS]
    require(set(values)=={0,1}and values.count(0)==84,'PR115 gap convention')
    points=[[1]+[indicator(p,k,j)for k,j in PAIRS]for p in PERMS]
    full_rank=rank(points);face_rank=rank([row for row,v in zip(points,values)if v==0])
    require((full_rank,face_rank)==(42,41),'exact dimension41 facet')
    gens=[]
    for j in range(4):
        swap=list(range(5));swap[j],swap[j+1]=swap[j+1],swap[j]
        gens.append([LOOKUP[tuple(swap[p[i]]for i in range(5))]for p in PERMS])
        gens.append([LOOKUP[tuple(p[swap[i]]for i in range(5))]for p in PERMS])
    gens.append([LOOKUP[tuple(p.index(i)for i in range(5))]for p in PERMS])
    initial=sum(1<<i for i,v in enumerate(values)if v)
    orbit={initial};queue=[initial]
    for mask in queue:
        support=[i for i in range(120)if mask>>i&1]
        for g in gens:
            transformed=sum(1<<g[i]for i in support)
            if transformed not in orbit:orbit.add(transformed);queue.append(transformed)
    require(len(orbit)==900,'full physical symmetry orbit')
    archive=json.loads((HERE.parent/'facets_n5.json').read_text())
    source=next(r for r in archive['orbits']if not r['minor']and r['vertices']==84)
    source_support=sum(1<<i for i in range(120)if i not in source['zero_set'])
    require(source_support in orbit,'same facet orbit as PR115')
    return {'vertices':120,'affine_dimension':full_rank-1,'face_dimension':face_rank-1,
            'zero_vertices':84,'positive_vertices':36,'physical_symmetry_orbit':len(orbit)}


def verify():
    results={name:gap_data(u)for name,u in fixtures().items()}
    require(results['block_zero_contact']['block_zero']and not results['block_zero_contact']['column_supported'],'exclusive block contact')
    require(results['column_supported_contact']['column_supported']and not results['column_supported_contact']['block_zero'],'exclusive column contact')
    require(results['zero_rhs_contact']['zero_rhs']and F(results['zero_rhs_contact']['tau'])>0,'exclusive zero-rhs contact')
    require(all(results[name]['connected_row_overlap']for name in results),'non-block-diagonal fixtures')
    require(F(results['strict_dense']['gap'])>0,'strict physical fixture')
    return {'facet':facet_data(),'fixtures':results,
            'theorem':'orbit A is nonnegative on every U(5); all900 symmetry copies inherit it',
            'limits':'not PMSH5, not a classification of rigid spectral8-contact sets'}


if __name__=='__main__':
    result=verify()
    if '--write-report'in sys.argv:(HERE/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    else:require(result==json.loads((HERE/'report.json').read_text()),'archived report')
    print(json.dumps(result,indent=2))
