#!/usr/bin/env python3
"""Replay six exact moving-fourth angular-neighborhood certificates."""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import hashlib,json
import neighborhood as n
c=n.c
SOURCE=c.HERE.parent/'complex-direction-cover'/'patches.json'
INDICES=(0,1,2,4,5,6)


def verify(certificate):
    source=SOURCE.read_bytes()
    if certificate['sourceSha256']!=hashlib.sha256(source).hexdigest():raise ValueError('Archived angular source changed')
    if [item['index']for item in certificate['patches']]!=list(INDICES):raise ValueError('Expected exactly the six quadratic sectors')
    patches=json.loads(source)
    for recorded in certificate['patches']:
        actual=n.certify(patches[recorded['index']]);actual['index']=recorded['index']
        if recorded!=actual:raise ValueError('Stored invariant-neighborhood bound changed')
    common=min(F(item['radius'])for item in certificate['patches'])
    if F(certificate['commonRadius'])!=common:raise ValueError('Common radius changed')
    return {'certifiedPatches':6,'strictlyPositiveTerms':384,'commonRadius':str(common),
            'invariantKernelDimensions':[p['invariantKernelDimension']for p in certificate['patches']],
            'scope':'Both angular charts k in [-1,1/2], a1 and a3 arbitrary complex; the two missing sectors are not asserted here.'}


def exact(value):
    if isinstance(value,bool)or not isinstance(value,(int,F)):raise TypeError('Use exact integer or Fraction coordinates')
    return F(value)


def coefficients_at(index,real,imag,certificate=None):
    real,imag=exact(real),exact(imag)
    certificate=json.loads((c.HERE/'certificate.json').read_text())if certificate is None else certificate
    bound=next((item for item in certificate['patches']if item['index']==index),None)
    if bound is None:raise ValueError('Uncertified patch index')
    if max(abs(real),abs(imag))>F(bound['radius']):raise ValueError('Outside the certified fourth-eigenvalue square')
    patch=json.loads(SOURCE.read_text())[index]
    keys,_,_=n.data(patch)
    columns=[c.column(patch,term,(real,imag))for term in patch['terms']]
    if any(set(col)-set(keys)for col in columns):raise ValueError('Unexpected coefficient row')
    matrix=[[col.get(key,F(0))for col in columns]for key in keys]
    rhs=[[F(key[0]=='mass')]for key in keys];selected=bound['selectedRows']
    values=c.solve_many([matrix[i]for i in selected],[rhs[i]for i in selected])
    if c.matmul(matrix,values)!=rhs:raise ValueError('Deformed coefficient identity failed')
    coefficients=[row[0]for row in values]
    if min(coefficients)<F(bound['minimumRetainedCoefficient'])or min(coefficients)<=0:raise ValueError('Deformed positivity bound failed')
    return patch,coefficients


def weights_at(index,real,imag,slope,radius):
    slope,radius=exact(slope),exact(radius)
    patch,coefficients=coefficients_at(index,real,imag)
    lo,hi=map(F,patch['interval'])
    if not lo<=slope<=hi:raise ValueError('Outside certified angular patch')
    degree=patch['degree'];alpha=(slope-lo)/(hi-lo);weights=[F(0)]*24
    for term,coefficient in zip(patch['terms'],coefficients):
        control,word,centers,_=term
        value=coefficient*comb(degree,control)*alpha**control*(1-alpha)**(degree-control)
        for center in centers:value*=(radius-F(center))**2
        permutation=c.a.PERMS.index(tuple(int(v)-1 for v in word));weights[permutation]+=value
    total=sum(weights)
    if total<=0 or min(weights)<0:raise ValueError('Radial normalizer failed')
    return [v/total for v in weights]


if __name__=='__main__':
    print(json.dumps(verify(json.loads((c.HERE/'certificate.json').read_text())),indent=2))
