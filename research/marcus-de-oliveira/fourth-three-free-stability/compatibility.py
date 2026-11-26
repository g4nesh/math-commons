"""Exact translated-spectrum coefficient matrices for the archived angular laws."""
from fractions import Fraction as F
from pathlib import Path
from functools import lru_cache
import importlib.util,json
from math import gcd
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('stability_cover',HERE.parent/'complex-direction-cover'/'check.py')
c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
a=c.atlas


def matmul(a,b):
    columns=list(zip(*b));integer_rows=[];integer_cols=[]
    for values,target in ((a,integer_rows),(columns,integer_cols)):
        for row in values:
            denominator=1
            for value in row:denominator=denominator*value.denominator//gcd(denominator,value.denominator)
            target.append(([value.numerator*(denominator//value.denominator)for value in row],denominator))
    return [[F(sum(x*y for x,y in zip(row,col)),dr*dc)for col,dc in integer_cols]for row,dr in integer_rows]


def solve_many(matrix,right):
    """Fraction-free elimination, integer back substitution, exact residual."""
    n=len(matrix)
    if not n or len(right)!=n or any(len(row)!=n for row in matrix):raise ValueError('Invalid exact system dimensions')
    width=len(right[0])
    if not width or any(len(row)!=width for row in right):raise ValueError('Invalid exact RHS dimensions')
    matrix=[list(map(F,row))for row in matrix];right=[list(map(F,row))for row in right]
    rows=[]
    for row,rhs in zip(matrix,right):
        values=row+rhs;denominator=1
        for value in values:denominator=denominator*value.denominator//gcd(denominator,value.denominator)
        rows.append([value.numerator*(denominator//value.denominator)for value in values])
    previous=1
    for k in range(n-1):
        pivot=next((i for i in range(k,n)if rows[i][k]),None)
        if pivot is None:raise ValueError('Singular exact system')
        rows[k],rows[pivot]=rows[pivot],rows[k];diagonal=rows[k][k]
        for i in range(k+1,n):
            factor=rows[i][k]
            for j in range(k+1,n+width):
                value=diagonal*rows[i][j]-factor*rows[k][j]
                quotient,remainder=divmod(value,previous)
                if remainder:raise ValueError('Nonexact fraction-free division')
                rows[i][j]=quotient
            rows[i][k]=0
        previous=diagonal
    determinant=rows[-1][n-1]
    if not determinant:raise ValueError('Singular exact system')
    # By Cramer's rule, determinant times every solution coordinate is an integer.
    scaled=[[0]*width for _ in range(n)]
    for i in range(n-1,-1,-1):
        for j in range(width):
            value=determinant*rows[i][n+j]-sum(rows[i][k]*scaled[k][j]for k in range(i+1,n))
            quotient,remainder=divmod(value,rows[i][i])
            if remainder:raise ValueError('Nonexact integer back substitution')
            scaled[i][j]=quotient
    answer=[[F(value,determinant)for value in row]for row in scaled]
    if matmul(matrix,answer)!=right:raise ValueError('Exact residual')
    return answer


@lru_cache(None)
def gaps(real=F(0),imag=F(0)):
    spectrum=tuple(a.add(b,(real,imag))for b in a.SPECTRUM)
    u=[[F(v,625)for v in row]for row in a.U_NUMERATORS]
    b=[[a.ZERO for _ in range(4)]for _ in range(4)]
    for i in range(4):
        for j in range(4):
            for k in range(4):b[i][j]=a.add(b[i][j],a.scale(u[i][k]*u[j][k],spectrum[k]))
    variables=[{(1,0,0):a.ONE},{(0,1,0):a.ONE},{(0,0,1):a.ONE},{}]
    matrix=[[a.poly_add(a.constant(b[i][j]),variables[i]if i==j else{})for j in range(4)]for i in range(4)]
    target=a.determinant(matrix)
    vertices=[a.poly_product(a.poly_add(variables[p[j]],a.constant(spectrum[j]))for j in range(4))for p in a.PERMS]
    result=tuple(a.poly_add(v,a.poly_scale(-1,target))for v in vertices)
    if any(poly.get((0,0,0),a.ZERO)!=a.ZERO for poly in result):raise ValueError('Scalar cancellation failed')
    return result


def column(patch,term,delta=(F(0),F(0))):
    chart=patch['chart'];lo,hi=map(F,patch['interval']);degree=patch['degree']
    control,word,origin,_=term;centers=tuple(map(F,origin));permutation=a.PERMS.index(tuple(int(v)-1 for v in word))
    entries={('mass',control):F(1)}
    for j,mult,slope in ((control,degree+1-control,lo),(control+1,control+1,hi)):
        omega=(F(1),slope)if chart=='h'else(slope,F(1))
        for (x,middle,y),value in gaps(*delta)[permutation].items():
            if middle:value=a.mul(value,omega)
            for radial,factor in enumerate(c.square_product(centers)):
                for coord in range(2):
                    key=(j,x,y,middle+radial,coord)
                    entries[key]=entries.get(key,F(0))+mult*factor*value[coord]
    return {key:v for key,v in entries.items()if v}


def matrices(patch):
    points=((F(0),F(0)),(F(1),F(0)),(F(-1),F(0)),(F(0),F(1)),(F(0),F(-1)))
    cols=[[column(patch,t,point)for t in patch['terms']]for point in points]
    keys=sorted(set().union(*(set(col)for group in cols for col in group)),key=repr)
    raw=[[[col.get(k,F(0))for col in group]for k in keys]for group in cols]
    base=raw[0];dx=[[(x-y)/2 for x,y in zip(r,s)]for r,s in zip(raw[1],raw[2])];dy=[[(x-y)/2 for x,y in zip(r,s)]for r,s in zip(raw[3],raw[4])]
    xx=[[(x+y)/2-z for x,y,z in zip(r,s,t)]for r,s,t in zip(raw[1],raw[2],base)]
    yy=[[(x+y)/2-z for x,y,z in zip(r,s,t)]for r,s,t in zip(raw[3],raw[4],base)]
    return keys,base,dx,dy,xx,yy


def independent_rows(rows):
    prime=1000000007;basis=[];indices=[]
    for i,row in enumerate(rows):
        reduced=[v.numerator*pow(v.denominator,-1,prime)%prime for v in row]
        for pivot,base in basis:
            factor=reduced[pivot];reduced=[(x-factor*y)%prime for x,y in zip(reduced,base)]
        pivot=next((j for j,v in enumerate(reduced)if v),None)
        if pivot is not None:
            inv=pow(reduced[pivot],-1,prime);basis.append((pivot,[v*inv%prime for v in reduced]));indices.append(i)
    return indices


def audit(patch):
    c.verify_patch(patch)
    keys,base,dx,dy,xx,yy=matrices(patch);weights=[F(t[3])for t in patch['terms']]
    rhs=[F(key[0]=='mass')for key in keys]
    assert all(sum(x*y for x,y in zip(row,weights))==b for row,b in zip(base,rhs))
    selected=independent_rows(base)
    if len(selected)!=len(weights):raise ValueError('Archived support not independent')
    derivative=[[sum(x*y for x,y in zip(row,weights))for row in matrix]for matrix in(dx,dy)]
    correction=solve_many([base[i]for i in selected],[[-derivative[j][i]for j in range(2)]for i in selected])
    residuals=[[sum(base[i][k]*correction[k][j]for k in range(len(weights)))+derivative[j][i]for i in range(len(keys))]for j in range(2)]
    return {'chart':patch['chart'],'interval':patch['interval'],'rows':len(keys),'columns':len(weights),'selectedRows':selected,
            'derivativeResiduals':[[{'row':list(keys[i]),'value':str(v)}for i,v in enumerate(values)if v]for values in residuals]},(keys,base,dx,dy,xx,yy,weights,rhs,selected,correction)


if __name__=='__main__':
    patches=json.loads((HERE.parent/'complex-direction-cover'/'patches.json').read_text())
    for index,patch in enumerate(patches):
        if len(patch['terms'])>64:continue
        result,_=audit(patch)
        print(json.dumps({'index':index,**result}),flush=True)
