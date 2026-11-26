#!/usr/bin/env python3
"""Exact complex rectangle with unrestricted first, third and fourth entries."""
from fractions import Fraction as F
from functools import lru_cache
from math import comb
from pathlib import Path
import argparse
import importlib.util
import json

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('prior_tube', HERE.parent/'two-middle-positive-laws/check.py')
g = importlib.util.module_from_spec(spec)
spec.loader.exec_module(g)
Q, Z, O, sp, PERMS = g.Q, g.Z, g.O, g.sp, g.PERMS
ALLOWED = tuple(j for j, p in enumerate(PERMS) if p[3] != 1)
MONOMIALS = tuple((a, b) for a in range(3) for b in range(3-a))
COLUMNS = tuple((p, e) for p in ALLOWED for e in MONOMIALS)
KEYS = tuple((r, e) for r, e in g.KEYS if r != 0)
LOWER, UPPER = F(1,4), F(3,4)


def load(name):
    return json.loads((HERE/name).read_text())


def norm(a):
    return max(sum(abs(v) for v in row) for row in a)


def coefficient_matrix(eta):
    output = []
    for row, (a,b) in KEYS:
        part, index = ('real',row-1) if row <= 4 else ('imag',row-5)
        free = g.prior.FREE[index]
        line = []
        for p, (c,d) in COLUMNS:
            gap = g.prior.gaps()[p]
            constant = gap[free]+eta*gap[free|{1}]
            slope = gap[free|{3}]+eta*gap[free|{1,3}]
            line.append(getattr(constant,part)*int((a,b)==(c,d))
                        +getattr(slope,part)*int((a,b)==(c+1,d))
                        +getattr(Q(0,1)*slope,part)*int((a,b)==(c,d+1)))
        output.append(line)
    return output


def center_controls(data):
    if data['degree'] != 3 or (F(data['lower']),F(data['upper'])) != (LOWER,UPPER):
        raise ValueError('Wrong center interval or degree')
    controls = [[F(0)]*108 for _ in range(4)]
    floors = {(k,p,v):F(0) for k in range(4) for p in ALLOWED
              for v in ((1,0,0),(0,1,0),(0,0,1))}
    lookup = {key:j for j,key in enumerate(COLUMNS)}
    seen = set()
    for k,p,vector,text in data['terms']:
        if type(k) is not int or k not in range(4) or type(p) is not int or p not in ALLOWED or len(vector)!=3 or any(type(v) is not int or v not in (-1,0,1) for v in vector) or type(text) is not str:
            raise ValueError('Malformed real-law term')
        key = (k,p,tuple(vector));mass=F(text)
        if not any(vector) or mass<=0 or key in seen:
            raise ValueError('Nonpositive or duplicate square term')
        seen.add(key)
        if key in floors:
            floors[key]+=mass
        for exponent,value in g.square(vector).items():
            controls[k][lookup[p,exponent]]+=mass*value
    if min(floors.values()) < 1:
        raise ArithmeticError('A positive-definite center floor is missing')
    return controls


def verify_relations(data,a0,ar,ai):
    exponents = [tuple(e) for e in data['exponents']]
    if len(set(exponents))!=len(exponents) or any(len(e)!=2 or any(type(k)is not int or k<0 for k in e) or sum(e)>3 for e in exponents):
        raise ValueError('Invalid relation exponents')
    if len(data['relations'])!=26:
        raise ValueError('Expected 26 polynomial relations')
    evaluated=[]
    for relation in data['relations']:
        residual={};at=[F(0)]*80;seen=set()
        for block,row,text in relation:
            if type(block)is not int or not 0<=block<len(exponents) or type(row)is not int or not 0<=row<80 or type(text)is not str or (block,row)in seen:
                raise ValueError('Malformed sparse relation')
            seen.add((block,row));value=F(text);u,v=exponents[block]
            if v==0:at[row]+=value*F(1,2)**u
            for du,dv,matrix in ((0,0,a0),(1,0,ar),(0,1,ai)):
                target=residual.setdefault((u+du,v+dv),[F(0)]*108)
                for j,x in enumerate(matrix[row]):
                    if x:target[j]+=value*x
        if any(x for row in residual.values() for x in row):
            raise ArithmeticError('A polynomial left-kernel identity fails')
        evaluated.append(at)
    columns=g.indices(data['relationColumns'],26,80)
    g.inverse([[row[j] for j in columns] for row in evaluated])
    # Independence at one point proves independence over Q(u,v), hence global rank <=54.
    return True


def verify(real_law,rank_data,bounds):
    if sp.e.multiply(sp.U,sp.e.adjoint(sp.U))!=sp.e.identity(4):
        raise ArithmeticError('Fixture is not unitary')
    controls=center_controls(real_law)
    a0=coefficient_matrix(Q(LOWER))
    ar=g.difference(coefficient_matrix(Q(UPPER)),a0)
    ai=g.difference(coefficient_matrix(Q(LOWER,1)),a0)
    end=[[x+y for x,y in zip(row,other)] for row,other in zip(a0,ar)]
    for j in range(5):
        left=g.matvec(a0,controls[j]) if j<4 else [F(0)]*80
        right=g.matvec(end,controls[j-1]) if j else [F(0)]*80
        if any((4-j)*x+j*y for x,y in zip(left,right)):
            raise ArithmeticError('Real-interval Bernstein identity fails')
    verify_relations(rank_data,a0,ar,ai)
    rows=g.indices(bounds['gapRows'],54,80);columns=g.indices(bounds['gapColumns'],54,108)
    select=lambda a:[[a[i][j]for j in columns]for i in rows]
    bc=[[x+F(1,2)*y for x,y in zip(row,other)]for row,other in zip(select(a0),select(ar))]
    inv=g.inverse(bc)
    x=g.matmul(inv,select(ar));y=g.matmul(inv,select(ai))
    forcing=g.matmul(inv,[[sum(ai[i][j]*controls[k][j]for j in range(108))for k in range(4)]for i in rows])
    power=[[F(i==j)for j in range(54)]for i in range(54)]
    row_sums=[F(0)]*54
    k=bounds['power']
    if type(k)is not int or not 1<=k<=64:
        raise ValueError('Invalid finite power')
    for j in range(k):
        for i,row in enumerate(power):row_sums[i]+=sum(abs(v)for v in row)/2**j
        power=g.matmul(power,x)
    tail=norm(power)/2**k
    if tail>=1:
        raise ArithmeticError('Real-interval power is not contractive')
    inverse_bound=max(row_sums)/(1-tail)
    imag_bound=norm(y)
    forcing_bound=max(abs(value)for row in forcing for value in row)
    if tail>=F(bounds['tailUpperBound']) or inverse_bound>=F(bounds['inverseBound']) or imag_bound>=F(bounds['imaginaryOperatorBound']) or forcing_bound>=F(bounds['forcingBound']):
        raise ArithmeticError('A declared rational upper bound fails')
    height=F(bounds['height']);K=F(bounds['inverseBound']);N=F(bounds['imaginaryOperatorBound']);M=F(bounds['forcingBound'])
    theta=height*K*N
    displacement=2*height*K*M/(1-theta)
    if height<=0 or theta>=F(1,2) or displacement>=F(4,5):
        raise ArithmeticError('Complex positivity estimate fails')
    energy=[F(0)]*4
    for k,p,vector,mass in real_law['terms']:energy[k]+=F(mass)*sum(v*v for v in vector)
    maximum=max(energy);energy_upper=(maximum.numerator+maximum.denominator-1)//maximum.denominator
    mass_floor=F(1,5*energy_upper+72)
    return dict(a0=a0,ar=ar,ai=ai,rows=rows,columns=columns,controls=controls,height=height,
                bc=bc,x=x,y=y,forcing=forcing,power=bounds['power'],tailUpperBound=F(bounds['tailUpperBound']),
                inverseBound=K,imaginaryOperatorBound=N,forcingBound=M,imaginaryContraction=theta,
                numeratorDisplacement=displacement,energyUpperBound=energy_upper,massFloor=mass_floor)



@lru_cache(None)
def prepared():
    return verify(load('real-law.json'),load('rank.json'),load('bounds.json'))


def weights_at(real_second,imag_second,real_fourth,imag_fourth):
    u,v,s,t=map(g.exact,(real_second,imag_second,real_fourth,imag_fourth))
    data=prepared()
    if not LOWER<=u<=UPPER or abs(v)>data['height']:
        raise ValueError('Second eigenvalue is outside the certified rectangle')
    tau=(u-LOWER)/(UPPER-LOWER)
    bern=[F(comb(3,k))*tau**k*(1-tau)**(3-k)for k in range(4)]
    coefficients=[sum(bern[k]*data['controls'][k][j]for k in range(4))for j in range(108)]
    matrix=[[data['a0'][i][j]+tau*data['ar'][i][j]+v*data['ai'][i][j]for j in data['columns']]for i in data['rows']]
    rhs=[[-v*sum(data['ai'][i][j]*coefficients[j]for j in range(108))]for i in data['rows']]
    correction=g.solve(matrix,rhs)
    for j,delta in zip(data['columns'],correction):coefficients[j]+=delta[0]
    numerators=[F(0)]*24
    for value,(p,(a,b))in zip(coefficients,COLUMNS):numerators[p]+=value*s**a*t**b
    denominator=sum(numerators)
    if denominator<=0 or any(numerators[p]<(1+s*s+t*t)/5 for p in ALLOWED):
        raise ArithmeticError('An exact evaluated numerator violates its proven floor')
    return [q/denominator for q in numerators]


def report():
    data=prepared();r=data['height'];checks=[]
    for u,v,s,t in ((LOWER,r,100,-100),(UPPER,-r,0,1),(F(1,2),r,F(7,3),F(-11,5)),(F(2,5),-r/3,-5,7)):
        weights=weights_at(u,v,s,t)
        value=g.prior.direct_check([Q(-7,3),Q(u,v),Q(11,-5),Q(s,t)],weights)
        checks.append({'second':[str(u),str(v)],'fourth':[str(s),str(t)],'determinant':value.encode()})
    return {'schemaVersion':1,'secondRealInterval':[str(LOWER),str(UPPER)],'secondImaginaryHeight':str(r),
            'unrestrictedComplexEntries':[1,3,4],'polynomialRelations':26,'globalRankUpperBound':54,
            'realLawTerms':len(load('real-law.json')['terms']),'power':data['power'],
            'bounds':{key:str(data[key])for key in ('tailUpperBound','inverseBound','imaginaryOperatorBound','forcingBound','imaginaryContraction','numeratorDisplacement')},
            'positiveNumeratorFloor':'(1+|a4|^2)/5','uniformMassLowerBound':str(data['massFloor']),'centerEnergyUpperBound':data['energyUpperBound'],'allowedPermutations':18,'zeroPermutations':6,'directChecks':checks}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-report',action='store_true');args=parser.parse_args()
    result=report();text=json.dumps(result,indent=2)+'\n'
    if args.write_report:(HERE/'report.json').write_text(text)
    elif load('report.json')!=result:raise ArithmeticError('Archived report differs from exact replay')
    print(text,end='')


if __name__=='__main__':main()
