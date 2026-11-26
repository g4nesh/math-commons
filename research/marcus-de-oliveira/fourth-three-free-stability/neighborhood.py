"""Exact invariant-subspace and positivity bounds for translated angular weights."""
from fractions import Fraction as F
from math import gcd
from pathlib import Path
import json
import compatibility as c


def lcm(a,b):return abs(a*b)//gcd(a,b)


def matmul(a,b):
    columns=list(zip(*b));rows=[];cols=[]
    for row in a:
        d=1
        for v in row:d=lcm(d,v.denominator)
        rows.append(([v.numerator*(d//v.denominator)for v in row],d))
    for col in columns:
        d=1
        for v in col:d=lcm(d,v.denominator)
        cols.append(([v.numerator*(d//v.denominator)for v in col],d))
    return [[F(sum(x*y for x,y in zip(row,col)),dr*dc)for col,dc in cols]for row,dr in rows]


def data(patch):
    keys,base,dx,dy,xx,yy=c.matrices(patch)
    columns=[c.column(patch,t,(F(1),F(1)))for t in patch['terms']]
    if any(set(col)-set(keys)for col in columns):raise ValueError('Missing mixed-coordinate row')
    xy=[[col.get(key,F(0))-base[i][j]-dx[i][j]-dy[i][j]-xx[i][j]-yy[i][j]for j,col in enumerate(columns)]for i,key in enumerate(keys)]
    if yy!=[[-v for v in row]for row in xx]:raise ValueError('Holomorphic quadratic coefficient mismatch')
    return keys,base,(dx,dy,xx,xy,yy)


def certify(patch):
    c.c.verify_patch(patch)
    keys,base,coefficients=data(patch);n=len(patch['terms']);weights=[F(term[3])for term in patch['terms']]
    rhs=[[F(key[0]=='mass')]for key in keys]
    if matmul(base,[[v]for v in weights])!=rhs:raise ValueError('Base identity changed')
    selected=c.independent_rows(base)
    if len(selected)!=n:raise ValueError('Base columns are not independent')
    combined=c.solve_many([base[i]for i in selected],[sum((matrix[i]for matrix in coefficients),[])for i in selected])
    operators=[[[row[k*n+j]for j in range(n)]for row in combined]for k in range(5)]
    other=[i for i in range(len(keys))if i not in set(selected)]
    residuals=[]
    for operator,matrix in zip(operators,coefficients):
        projected=matmul([base[i]for i in other],operator)
        residuals.extend([[x-y for x,y in zip(row,matrix[i])]for row,i in zip(projected,other)])
    independent=c.independent_rows(residuals);basis=[residuals[i]for i in independent]
    if not basis:raise ValueError('Expected a nontrivial residual space')
    columns=c.independent_rows([list(row)for row in zip(*basis)])
    if len(columns)!=len(basis):raise ValueError('Residual basis rank mismatch')
    # Normalize pivot columns to identity: H = R_pivot^-1 R.
    pivot=[[row[j]for j in columns]for row in basis]
    normalized=c.solve_many(pivot,basis)
    if [[row[j]for j in columns]for row in normalized]!=[[F(i==j)for j in range(len(basis))]for i in range(len(basis))]:raise ValueError('Residual pivot normalization failed')
    if matmul([[row[j]for j in columns]for row in residuals],normalized)!=residuals:raise ValueError('Residual row span not exact')
    if any(v[0] for v in matmul(normalized,[[v]for v in weights])):raise ValueError('Base weights do not lie in invariant kernel')
    for operator in operators:
        transformed=matmul(normalized,operator)
        if matmul([[row[j]for j in columns]for row in transformed],normalized)!=transformed:raise ValueError('Residual kernel is not invariant')
    norms=[];actions=[]
    for operator in operators:
        norms.append(max(sum(abs(v)*weights[j]/weights[i]for j,v in enumerate(row))for i,row in enumerate(operator)))
        actions.append(max(abs(v[0])/weights[i]for i,v in enumerate(matmul(operator,[[w]for w in weights]))))
    radius=F(1)
    while True:
        rho=radius*sum(norms[:2])+radius*radius*sum(norms[2:])
        forcing=radius*sum(actions[:2])+radius*radius*sum(actions[2:])
        if rho<1 and 2*forcing<=1-rho:break
        radius/=2
    result={'chart':patch['chart'],'interval':patch['interval'],'positiveTerms':n,'coefficientRows':len(keys),
            'residualRank':len(basis),'invariantKernelDimension':n-len(basis),'selectedRows':selected,
            'radius':str(radius),'operatorNorms':list(map(str,norms)),'operatorActions':list(map(str,actions)),
            'rho':str(rho),'forcing':str(forcing),'relativeChangeBound':str(forcing/(1-rho)),
            'minimumRetainedCoefficient':str(min(weights)/2)}
    return result


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--index',type=int,default=0);parser.add_argument('--all',action='store_true');parser.add_argument('--report',type=Path);args=parser.parse_args()
    patches=json.loads((c.HERE.parent/'complex-direction-cover/patches.json').read_text())
    indices=[i for i,p in enumerate(patches)if len(p['terms'])==64]if args.all else[args.index]
    report=[]
    for i in indices:
        result=certify(patches[i]);result['index']=i;report.append(result)
        print(json.dumps({'index':i,'radius':result['radius'],'residualRank':result['residualRank']}),flush=True)
    if args.report:args.report.write_text(json.dumps(report,indent=2)+'\n')
