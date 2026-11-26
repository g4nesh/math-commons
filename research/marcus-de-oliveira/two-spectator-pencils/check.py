#!/usr/bin/env python3
"""Exact physical two-variable certificate; only Python standard library."""
from fractions import Fraction as F
from itertools import permutations
import json

# (x,y) denotes x+y*j with j^2=-3.
def q(x=0,y=0):return F(x),F(y)
ZERO,ONE=q(),q(1)
def add(x,y):return x[0]+y[0],x[1]+y[1]
def mul(x,y):return x[0]*y[0]-3*x[1]*y[1],x[0]*y[1]+x[1]*y[0]
def scale(t,x):return t*x[0],t*x[1]
def conj(x):return x[0],-x[1]
def poly_add(a,b):
    out=dict(a)
    for key,value in b.items():out[key]=add(out.get(key,ZERO),value)
    return {key:value for key,value in out.items() if value!=ZERO}
def poly_mul(a,b):
    out={}
    for (i,j),x in a.items():
        for (r,s),y in b.items():out[i+r,j+s]=add(out.get((i+r,j+s),ZERO),mul(x,y))
    return {key:value for key,value in out.items() if value!=ZERO}
def constant(x):return {} if x==ZERO else {(0,0):x}
def poly_det(matrix):
    total={};n=len(matrix)
    for p in permutations(range(n)):
        term={(0,0):ONE}
        for i,j in enumerate(p):term=poly_mul(term,matrix[i][j])
        sign=(-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
        total=poly_add(total,{key:scale(sign,value) for key,value in term.items()})
    return total
def adjoint(a):return [[conj(a[j][i]) for j in range(len(a))] for i in range(len(a))]
def matrix_mul(a,b):
    return [[sum_q(mul(a[i][k],b[k][j]) for k in range(len(b))) for j in range(len(b[0]))] for i in range(len(a))]
def sum_q(xs):
    result=ZERO
    for x in xs:result=add(result,x)
    return result

def main():
    G=[[q(F(x,2)) for x in row] for row in [[1,-1,1],[1,1,-1],[-1,1,1]]]
    assert matrix_mul(G,adjoint(G))==matrix_mul(adjoint(G),G)
    spectrum=[q(F(1,2)),q(F(1,2),F(1,2)),q(F(1,2),F(-1,2))]
    for eigenvalue in spectrum:
        characteristic=[[constant(add(x,scale(-1,eigenvalue)) if i==j else x) for j,x in enumerate(row)] for i,row in enumerate(G)]
        assert poly_det(characteristic)=={}
    matrix=[[{} for _ in range(4)] for _ in range(4)]
    for i in range(3):
        for j in range(3):
            matrix[i][j]=poly_add(matrix[i][j],constant(G[i][j]))
            matrix[i+1][j+1]=poly_add(matrix[i+1][j+1],constant(G[i][j]))
    matrix[3][3]=poly_add(matrix[3][3],{(1,0):ONE})
    matrix[0][0]=poly_add(matrix[0][0],{(0,1):ONE})
    actual=poly_det(matrix)
    expected={(0,0):q(1),(1,0):q(F(3,2)),(0,1):q(F(3,2)),(1,1):q(2)}
    assert actual==expected
    A=[{(1,0):ONE}]+[constant(x) for x in spectrum]
    B=[{(0,1):ONE}]+[constant(x) for x in spectrum]
    weights=[((0,2,3,1),F(2,7)),((1,0,2,3),F(1,3)),((2,3,1,0),F(4,21)),((3,2,0,1),F(4,21))]
    assert sum(w for p,w in weights)==1 and all(w>0 for p,w in weights)
    mixed={}
    for p,w in weights:
        vertex={(0,0):ONE}
        for i,j in enumerate(p):vertex=poly_mul(vertex,poly_add(A[i],B[j]))
        mixed=poly_add(mixed,{key:scale(w,value) for key,value in vertex.items()})
    assert mixed==actual
    # h=-1, k=0: strict failure of the ellipse proxy, while the polynomial certificate holds.
    h=F(-1);outside=h+F(1,2);d=F(1,2)+1/(4*outside);e=F(1,2)-1/(4*outside)
    assert d==0 and e==1 and abs(e)>abs(d)
    assert F(1)+F(3,2)*h==F(-1,2)
    print(json.dumps({'status':'passed','field':'Q(sqrt(-3))','original_dimension':4,'physical_normal_blocks':True,'arbitrary_complex_parameters':['h','k'],'determinant_coefficients':{'1':'1','h':'3/2','k':'3/2','hk':'2'},'weights':[{'permutation_zero_based':list(p),'weight':str(w)} for p,w in weights],'proxy_failure_example':{'h':'-1','k':'0','determinant':'-1/2','proxy_target_modulus':'1/2','proxy_disk_radius':'1/4'},'evidence':'exact polynomial identity and positive probability weights; not general MOC'},indent=2))
if __name__=='__main__':main()
