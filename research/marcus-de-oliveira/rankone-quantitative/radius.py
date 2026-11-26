#!/usr/bin/env python3
"""Exact quantitative moving-fourth radius for the two-free-eigenvalue identity."""
from fractions import Fraction as F
from pathlib import Path
import importlib.util,sys,json
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('quant_rankone',HERE.parent/'two-free-rank-one'/'check.py')
c=importlib.util.module_from_spec(spec);sys.modules[spec.name]=c;spec.loader.exec_module(c)
SOURCE=json.loads((HERE.parent/'two-free-rank-one'/'certificate.json').read_text())
IDS={''.join(str(i+1)for i in p):j for j,p in enumerate(c.PERMS)}
PAIRS=tuple((IDS[a],IDS[b])for a,b,_ in SOURCE['products'])
BASE=tuple(F(v)for a,b,v in SOURCE['products'])+tuple(map(F,SOURCE['phase']))
RHS=[F(0)]*36+[F(1)]


def matmul(a,b):
    columns=list(zip(*b))
    return [[sum(x*y for x,y in zip(row,col))for col in columns]for row in a]

def solve_many(a,b):
    n=len(a);width=len(b[0]);r=[list(map(F,row))+list(map(F,right))for row,right in zip(a,b)]
    if any(len(row)!=n for row in a)or len(b)!=n or any(len(row)!=width for row in b):raise ValueError('Invalid square solve dimensions')
    for j in range(n):
        k=next((i for i in range(j,n)if r[i][j]),None)
        if k is None:raise ValueError('Singular certificate matrix')
        r[j],r[k]=r[k],r[j];pivot=r[j][j];r[j]=[x/pivot for x in r[j]]
        for i in range(n):
            if i!=j and r[i][j]:
                factor=r[i][j];r[i]=[x-factor*y for x,y in zip(r[i],r[j])]
    result=[row[n:]for row in r]
    if matmul(a,result)!=[[F(v)for v in row]for row in b]:raise ValueError('Exact linear solve residual')
    return result


def gap_data():
    base=c.gap_rows();rx=c.gap_rows(fourth=c.s.Q(1));ry=c.gap_rows(fourth=c.s.Q(0,1))
    return base,[[b-a for a,b in zip(r,s)]for r,s in zip(base,rx)],[[b-a for a,b in zip(r,s)]for r,s in zip(base,ry)]

def coefficient_matrices():
    rows,dx,dy=gap_data();columns=[[]for _ in range(6)]
    for i,j in PAIRS:
        p=c.quadratic_product
        pieces=[p(rows[i],rows[j]),[a+b for a,b in zip(p(dx[i],rows[j]),p(rows[i],dx[j]))],
                [a+b for a,b in zip(p(dy[i],rows[j]),p(rows[i],dy[j]))],p(dx[i],dx[j]),
                [a+b for a,b in zip(p(dx[i],dy[j]),p(dy[i],dx[j]))],p(dy[i],dy[j])]
        for k,col in enumerate(pieces):columns[k].append(col+[F(k==0)])
    real,imag=c.determinant_parts()
    for k,cols in enumerate(columns):cols.extend([[-F(v)if k==0 else F(0)for v in real]+[F(0)],[-F(v)if k==0 else F(0)for v in imag]+[F(0)]])
    return [[[F(v)for v in row]for row in zip(*cols)]for cols in columns]

def matrix_at(x,y,matrices=None):
    matrices=coefficient_matrices()if matrices is None else matrices;coefficients=(1,x,y,x*x,x*y,y*y)
    return [[sum(a*matrix[i][j]for a,matrix in zip(coefficients,matrices))for j in range(37)]for i in range(37)]

def direct_matrix(x,y):
    rows=c.gap_rows(fourth=c.s.Q(x,y));columns=[c.quadratic_product(rows[i],rows[j])+[F(1)]for i,j in PAIRS]
    real,imag=c.determinant_parts();columns.extend([[-F(v)for v in real]+[F(0)],[-F(v)for v in imag]+[F(0)]])
    return [list(row)for row in zip(*columns)]

def derive_bound():
    c.verify(SOURCE)
    if len(PAIRS)!=35 or min(BASE[:35])<=0:raise ValueError('Expected 35 strictly positive base products')
    matrices=coefficient_matrices();a0=matrices[0]
    if matmul(a0,[[v]for v in BASE])!=[[v]for v in RHS]:raise ValueError('Base coefficient identity failed')
    identity=[[F(i==j)for j in range(37)]for i in range(37)]
    inverse=solve_many(a0,identity)
    if matmul(inverse,a0)!=identity:raise ValueError('Inverse residual failed')
    scales=list(BASE[:35])+[max(F(1),abs(v))for v in BASE[35:]]
    norms=[];actions=[]
    for matrix in matrices[1:]:
        preconditioned=matmul(inverse,matrix)
        norms.append(max(sum(abs(v)*scales[j]/scales[i]for j,v in enumerate(row))for i,row in enumerate(preconditioned)))
        actions.append(max(abs(sum(v*BASE[j]for j,v in enumerate(row)))/scales[i]for i,row in enumerate(preconditioned)))
    linear=norms[0]+norms[1];quadratic=sum(norms[2:]);action_linear=sum(actions[:2]);action_quadratic=sum(actions[2:])
    rows,dx,dy=gap_data();witness=list(map(F,SOURCE['strictSeparator']))
    margins=[sum(a*b for a,b in zip(witness,row))for row in rows]
    slopes=[abs(sum(a*b for a,b in zip(witness,rowx)))+abs(sum(a*b for a,b in zip(witness,rowy)))for rowx,rowy in zip(dx,dy)]
    margin=min(margins);slope=max(slopes)
    if margin!=F(SOURCE['minimumGap'])or margin<=0:raise ValueError('Strict separation margin changed')
    radius=F(1);exponent=0
    while True:
        rho=radius*linear+radius*radius*quadratic
        forcing=radius*action_linear+radius*radius*action_quadratic
        affine_margin=min(v-radius*d for v,d in zip(margins,slopes))
        if rho<1 and 2*forcing<=1-rho and affine_margin>=margin/2:break
        radius/=2;exponent+=1
        if exponent>256:raise ValueError('Radius discovery exceeded the finite bound')
    return {'baseSourceCommit':'cb253d9c3c7746fd2278db9d41a5eb1668ed52b4',
            'fixedThird':['2','-1/4'],'fourthRealInterval':[str(-radius),str(radius)],'fourthImaginaryInterval':[str(-radius),str(radius)],
            'dyadicRadiusExponent':exponent,'positiveProducts':35,'systemSize':37,
            'preconditionedCoefficientNorms':dict(zip(('x','y','xx','xy','yy'),map(str,norms))),
            'linearNormBound':str(linear),'quadraticNormBound':str(quadratic),'perturbationNormBound':str(rho),
            'actionCoefficientNorms':dict(zip(('x','y','xx','xy','yy'),map(str,actions))),
            'actionNormBound':str(forcing),'relativeSolutionChangeBound':str(forcing/(1-rho)),'minimumBaseProduct':str(min(BASE[:35])),
            'minimumProductThroughoutSquare':str(min(BASE[:35])/2),'baseUniformLawMargin':str(margin),
            'maximumSeparatorSlope':str(slope),'minimumAffineSeparatorMargin':str(affine_margin),'uniformLawMarginThroughoutSquare':str(margin/2),
            'scope':'Pointwise MDO for all complex first two eigenvalues, throughout this explicit complex fourth-eigenvalue square; no common two-variable law anywhere in the square.'}


def exact_parameters(x,y,radius):
    if any(isinstance(v,bool)or not isinstance(v,(int,F))for v in(x,y)):raise TypeError('Use exact integer or Fraction coordinates')
    x,y=F(x),F(y)
    if abs(x)>radius or abs(y)>radius:raise ValueError('Outside the certified fourth-eigenvalue square')
    return x,y

def coefficients_at(x,y,bound=None):
    bound=json.loads((HERE/'bound.json').read_text())if bound is None else bound
    radius=F(bound['fourthRealInterval'][1]);x,y=exact_parameters(x,y,radius)
    values=[row[0]for row in solve_many(direct_matrix(x,y),[[v]for v in RHS])]
    if min(values[:35])<=0 or min(values[:35])<F(bound['minimumProductThroughoutSquare'])or sum(values[:35])!=1:raise ValueError('Positive product certificate failed')
    rows=c.gap_rows(fourth=c.s.Q(x,y));witness=list(map(F,SOURCE['strictSeparator']))
    minimum_margin=min(sum(a*b for a,b in zip(witness,row))for row in rows)
    if minimum_margin<=0 or minimum_margin<F(bound['uniformLawMarginThroughoutSquare']):raise ValueError('Uniform-law separation bound failed')
    return values


def verify(bound):
    if derive_bound()!=bound:raise ValueError('Recorded quantitative bound differs from exact derivation')
    matrices=coefficient_matrices()
    for x,y in ((F(0),F(0)),(F(1),F(0)),(F(0),F(1)),(F(1),F(1)),(F(-2,3),F(5,7))):
        if matrix_at(x,y,matrices)!=direct_matrix(x,y):raise ValueError('Polynomial matrix construction failed')
    radius=F(bound['fourthRealInterval'][1]);checks=[]
    for x,y in ((F(0),F(0)),(-radius,-radius),(-radius,radius),(radius,-radius),(radius,radius)):
        values=coefficients_at(x,y,bound)
        checks.append({'fourth':[str(x),str(y)],'minimumProduct':str(min(values[:35])),'phaseNonzero':any(values[35:])})
    if c.s.e.multiply(c.s.U,c.s.e.adjoint(c.s.U))!=c.s.e.identity(4):raise ValueError('Archived U is not unitary')
    return {'bound':bound,'exactMatrixChecks':5,'cornerCertificates':checks,
            'proof':'Relative Neumann-series bound and affine separator margins; the finite corner checks are regression evidence, not the neighborhood proof.'}


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--derive',action='store_true');parser.add_argument('--report',type=Path)
    args=parser.parse_args();path=HERE/'bound.json'
    if args.derive:path.write_text(json.dumps(derive_bound(),indent=2)+'\n')
    report=verify(json.loads(path.read_text()))
    if args.report:args.report.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
