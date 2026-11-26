#!/usr/bin/env python3
"""Reconstruct the archived rational continuations from support alone, stdlib only."""
import json
from fractions import Fraction as F
import patches as p

def mod(x,prime=1000000007):return x.numerator*pow(x.denominator,-1,prime)%prime


def independent_rows(rows):
    prime=1000000007;basis=[];indices=[]
    for i,row in enumerate(rows):
        reduced=list(map(mod,row))
        for pivot,base in basis:
            factor=reduced[pivot];reduced=[(a-factor*b)%prime for a,b in zip(reduced,base)]
        pivot=next((j for j,v in enumerate(reduced)if v),None)
        if pivot is not None:
            inverse=pow(reduced[pivot],-1,prime)
            basis.append((pivot,[v*inverse%prime for v in reduced]));indices.append(i)
    return indices


def krylov_degree(matrix,vector):
    prime=1000000007;matrix=[list(map(mod,row))for row in matrix];vector=list(map(mod,vector));basis=[]
    for degree in range(len(vector)+1):
        reduced=vector[:]
        for pivot,row in basis:
            factor=reduced[pivot];reduced=[(a-factor*b)%prime for a,b in zip(reduced,row)]
        pivot=next((j for j,v in enumerate(reduced)if v),None)
        if pivot is None:return degree
        inverse=pow(reduced[pivot],-1,prime);basis.append((pivot,[v*inverse%prime for v in reduced]))
        vector=[sum(a*b for a,b in zip(row,vector))%prime for row in matrix]
    raise ValueError('Krylov degree exceeds dimension')


def reconstruct(patch,data=None):
    terms=[{'centers':t['centers'],'permutation':t['permutation']}for t in patch['terms']]
    a0,a1=p.matrices(terms,patch['chart'],F(patch['center']),data)
    rows=independent_rows(a0)
    if len(rows)!=len(terms):raise ValueError('Support not independent')
    solved=p.solve_many([a0[i]for i in rows],[[F(i==48)]+a1[i]for i in rows])
    weights=[row[0]for row in solved];change=[row[1:]for row in solved]
    if min(weights)<=0 or p.matmul(a0,[[v]for v in weights])!=[[F(i==48)]for i in range(49)]:raise ValueError('Invalid exact center weights')
    if p.matmul(a0,change)!=a1:raise ValueError('Angular column space changes')
    degree=krylov_degree(change,weights);sequence=[weights]
    for _ in range(degree):sequence.append([sum(x*y for x,y in zip(row,sequence[-1]))for row in change])
    krylov=[list(row)for row in zip(*sequence[:-1])];indices=independent_rows(krylov)
    if len(indices)!=degree:raise ValueError('Unexpected exact Krylov rank')
    relation=[row[0]for row in p.solve_many([krylov[i]for i in indices],[[-sequence[-1][i]]for i in indices])]
    if any(sum(a*b for a,b in zip(row,relation))+last!=0 for row,last in zip(krylov,sequence[-1])):raise ValueError('Krylov relation residual')
    denominator=[F(1)]+[(-1)**j*relation[degree-j]for j in range(1,degree+1)]
    numerators=[weights]
    for power in range(1,degree):
        numerators.append([denominator[power]*v-sum(a*b for a,b in zip(row,numerators[-1]))for v,row in zip(weights,change)])
    out=dict(patch,denominator=list(map(str,denominator)),terms=[dict(term,numerator=list(map(str,poly)))for term,poly in zip(terms,zip(*numerators))])
    if out!=patch:raise ValueError('Archived continuation differs from exact reconstruction')
    p.check_patch(out,data)
    return {'chart':patch['chart'],'interval':patch['interval'],'columns':len(terms),'krylovDegree':degree,'exactReconstructionMatches':True}


if __name__=='__main__':
    certificate=json.loads((p.HERE/'certificate.json').read_text());data=p.pencil_data()
    print(json.dumps([reconstruct(patch,data)for patch in certificate['patches']],indent=2))
