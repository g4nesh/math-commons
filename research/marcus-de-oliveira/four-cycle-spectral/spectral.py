#!/usr/bin/env python3
"""Exact interval and complex-plane certificates for one fixed four-cycle.

Standard library only. Polynomial coefficient identities prove continuous
claims; finite sample evaluations are extra regressions, not the proof.
"""
from fractions import Fraction as F
from itertools import permutations
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import resource
import sys

ROOT = Path(__file__).resolve().parents[1]
VERIFY = ROOT / 'certificates' / 'verify.py'
spec = importlib.util.spec_from_file_location('spectral_exact', VERIFY)
e = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = e
spec.loader.exec_module(e)
Q, Z, O, I = e.QComplex, e.ZERO, e.ONE, e.QComplex(0,1)
U = [[Q(F(v,625)) for v in row] for row in
     ((-31,492,-240,300),(-300,225,500,0),(240,-180,225,500),(-492,-256,-180,225))]
BASE_A = [Q(-2,-1),Q(-1),Q(0,1),Q(1,-1)]
B_SPECTRUM = [Q(3,-2),Q(2),Q(1,2),Q(0,-1)]
TRIANGLE = ((1,3,2,0),(0,2,1,3),(1,2,0,3))
PLANE = ((0,1,2,3),(0,1,3,2),(0,2,1,3),(1,0,2,3),(3,1,2,0))
INTERVAL = (F(-40),F(1))
NUMERATORS = ((13203996,-5027220,-131806),(14562000,-10814944,6202444),(19109004,-15407836,1741862))
DENOMINATOR = (46875000,-31250000,7812500)


def product(values):
    out = O
    for value in values:
        out *= value
    return out


def leibniz(matrix):
    out = Z
    for p in permutations(range(len(matrix))):
        sign = (-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))
        out += sign*product(matrix[i][p[i]] for i in range(len(p)))
    return out


def trim(poly):
    poly = list(poly)
    while len(poly)>1 and poly[-1] == 0:
        poly.pop()
    return tuple(poly)


def add(p,q):
    n = max(len(p),len(q))
    return tuple((p[i] if i<len(p) else 0)+(q[i] if i<len(q) else 0) for i in range(n))


def multiply(p,q):
    out = [0]*(len(p)+len(q)-1)
    for i,x in enumerate(p):
        for j,y in enumerate(q):
            out[i+j] += x*y
    return tuple(out)


def evaluate(poly,t):
    if isinstance(t,bool) or not isinstance(t,(int,F,Q)):
        raise TypeError('Exact integer, Fraction, or Gaussian rational required')
    out = poly[-1]
    for coefficient in reversed(poly[:-1]):
        out = out*t+coefficient
    return out


def subtract(p,q):
    return add(p,tuple(-x for x in q))


def cross(x,y):
    return x.real*y.imag-x.imag*y.real


def affine_cross(x,y):
    return (cross(x[0],y[0]),cross(x[0],y[1])+cross(x[1],y[0]),cross(x[1],y[1]))


def vertex(a,p):
    return product(a[p[j]]+B_SPECTRUM[j] for j in range(4))


def matrix_sum(a):
    orbit = e.multiply(e.multiply(U,e.diagonal(B_SPECTRUM)),e.adjoint(U))
    return [[orbit[i][j]+(a[i] if i==j else Z) for j in range(4)] for i in range(4)]


def affine_data(a0,direction,ps):
    a = BASE_A.copy(); a[0] = a0
    matrix = matrix_sum(a)
    constant = e.determinant(matrix)
    assert constant == leibniz(matrix)
    cofactor = e.determinant([row[1:] for row in matrix[1:]])
    assert cofactor == leibniz([row[1:] for row in matrix[1:]])
    target = (constant,direction*cofactor)
    points = []
    for p in ps:
        varying_column = p.index(0)
        points.append((vertex(a,p),direction*product(a[p[j]]+B_SPECTRUM[j] for j in range(4) if j != varying_column)))
    return target,tuple(points)


def exact_solve(a,b):
    n = len(b)
    rows = [list(map(F,row))+[F(value)] for row,value in zip(a,b)]
    for j in range(n):
        pivot = next((i for i in range(j,n) if rows[i][j]),None)
        if pivot is None:
            raise ValueError('Coefficient system is singular')
        rows[j],rows[pivot] = rows[pivot],rows[j]
        factor = rows[j][j]; rows[j] = [x/factor for x in rows[j]]
        for i in range(n):
            if i != j:
                factor = rows[i][j]
                rows[i] = [x-factor*y for x,y in zip(rows[i],rows[j])]
    return tuple(row[-1] for row in rows)


def interval_algebra():
    target,points = affine_data(BASE_A[0],I,TRIANGLE)
    p,q,r = points
    area = affine_cross(subtract(q,p),subtract(r,p))
    raw = (affine_cross(subtract(q,target),subtract(r,target)),
           affine_cross(subtract(r,target),subtract(p,target)),
           affine_cross(subtract(p,target),subtract(q,target)))
    numerators = tuple(tuple(F(x)*15625 for x in poly) for poly in raw)
    denominator = tuple(F(x)*15625 for x in area)
    assert numerators == NUMERATORS and denominator == DENOMINATOR
    assert add(add(numerators[0],numerators[1]),numerators[2]) == denominator
    lhs = (Z,)
    for numerator,point in zip(numerators,points):
        lhs = add(lhs,multiply(numerator,point))
    assert lhs == multiply(denominator,target)
    return target,points,area,numerators,denominator


def positivity():
    left,right = INTERVAL
    n1,n2,n3 = NUMERATORS
    assert n1[2] < 0 and min(evaluate(n1,left),evaluate(n1,right)) > 0
    discriminant = n2[1]**2-4*n2[0]*n2[2]
    assert n2[2]>0 and discriminant<0
    assert n3[2]>0 and n3[1]+2*n3[2]*right<0 and evaluate(n3,right)>0
    # Denominator = 7812500*((t-2)^2+2), positive on all real t.
    assert DENOMINATOR == (6*7812500,-4*7812500,7812500)
    # Exact rational isolation of the first sign changes around t=0.
    assert evaluate(n1,F(-41))<0<evaluate(n1,F(-40))
    assert evaluate(n1,F(2))>0>evaluate(n1,F(3))
    assert evaluate(n3,F(1))>0>evaluate(n3,F(3,2))
    assert evaluate(n3,F(7))<0<evaluate(n3,F(8))
    return {'n1ConcaveEndpointValues':[str(evaluate(n1,t)) for t in INTERVAL],
            'n2Discriminant':str(discriminant),
            'n3DerivativeUpperBound':str(n3[1]+2*n3[2]*right),
            'n3RightEndpointValue':str(evaluate(n3,right)),
            'denominatorGlobalLowerBound':'15625000',
            'rootIsolation':{'N1NegativeRoot':['-41','-40'],'N1PositiveRoot':['2','3'],
                             'N3SmallerRoot':['1','3/2'],'N3LargerRoot':['7','8']}}


def interval_weights(t):
    if isinstance(t,bool) or not isinstance(t,(int,F)):
        raise TypeError('Interval parameters require exact integers or Fractions')
    if not INTERVAL[0]<=t<=INTERVAL[1]:
        raise ValueError('Parameter is outside the certified closed interval [-40,1]')
    return tuple(F(evaluate(p,t),evaluate(DENOMINATOR,t)) for p in NUMERATORS)


def plane_algebra():
    target,points = affine_data(Z,O,PLANE)
    rows = [[F(1)]*5]+[[p[k].real if part=='real' else p[k].imag for p in points]
                            for k in range(2) for part in ('real','imag')]
    rhs = [F(1),target[0].real,target[0].imag,target[1].real,target[1].imag]
    weights = exact_solve(rows,rhs)
    expected = (F(22689,390625),F(8848,15625),F(1296,15625),F(144,15625),F(110736,390625))
    assert weights == expected and sum(weights)==1 and min(weights)>0
    assert tuple(sum((w*p[k] for w,p in zip(weights,points)),Z) for k in range(2)) == target
    return target,points,weights


def check_at(a0,ps,weights):
    a = BASE_A.copy();a[0]=a0
    matrix = matrix_sum(a)
    actual = e.determinant(matrix)
    assert actual == leibniz(matrix)
    assert sum(weights)==1 and min(weights)>=0
    value = sum((w*vertex(a,p) for w,p in zip(weights,ps)),Z)
    assert value == actual
    return actual


def source_check():
    source = ROOT/'single-cycle'/'results.json'
    previous = json.loads(source.read_text())['fourCycleGreedyRuleFailure']
    assert [[z.encode() for z in row] for row in U] == previous['U']
    assert [z.encode() for z in BASE_A] == previous['determinant']['a']
    assert [z.encode() for z in B_SPECTRUM] == previous['determinant']['b']
    assert [list(p) for p in TRIANGLE] == [w['columnToRow'] for w in previous['determinant']['originalHullConvexWeights']]
    return source


def report():
    source = source_check()
    assert e.multiply(U,e.adjoint(U)) == e.identity(4) == e.multiply(e.adjoint(U),U)
    target,points,area,numerators,denominator = interval_algebra()
    sign_proof = positivity()
    checkpoints = []
    for t in (F(-40),F(-20),F(-1),F(0),F(1,2),F(1)):
        weights = interval_weights(t)
        value = check_at(BASE_A[0]+I*t,TRIANGLE,weights)
        checkpoints.append({'t':str(t),'weights':list(map(str,weights)),'determinant':value.encode()})
    assert interval_weights(0) != interval_weights(1)
    global_target,global_points,global_weights = plane_algebra()
    for z in (Z,O,Q(0,1),Q(-2,-1),Q(F(7,3),F(-11,5)),Q(100,-100)):
        check_at(z,PLANE,global_weights)
    encode_affine = lambda p:[z.encode() for z in p]
    return {'schemaVersion':1,'evidenceLevel':'written continuous parameter proofs and exact coefficient checks; same-owner, not Lean or independent review',
            'pythonVersion':platform.python_version(),'U':[[z.encode() for z in row]for row in U],
            'baseA':[z.encode() for z in BASE_A],'fixedB':[z.encode() for z in B_SPECTRUM],
            'deformation':'a0(t)=-2+i(t-1); every other eigenvalue and U fixed; t real',
            'interval':list(map(str,INTERVAL)),'supportColumnToRow':[list(p)for p in TRIANGLE],
            'determinantAffineConstantSlope':encode_affine(target),
            'verticesAffineConstantSlope':[encode_affine(p)for p in points],
            'orientedAreaCoefficients':list(map(str,area)),
            'weightNumeratorsAscending':[list(map(str,p)) for p in numerators],
            'weightCommonDenominatorAscending':list(map(str,denominator)),
            'polynomialIdentitiesVerified':True,'positivityProof':sign_proof,'exactParameterChecks':checkpoints,
            'complexPlaneExtension':{'claim':'For every complex a0=z, the other seven eigenvalues and U fixed, these five constant nonnegative weights reconstruct the determinant.',
                'supportColumnToRow':[list(p)for p in PLANE],'weights':list(map(str,global_weights)),
                'determinantAffineConstantSlope':encode_affine(global_target),
                'verticesAffineConstantSlope':[encode_affine(p)for p in global_points],
                'coefficientIdentitiesVerified':True,'discoveryCredit':'ganesh-coordinator; independently solved and checked here'},
            'scope':'One fixed four-cycle and seven fixed eigenvalues. No general MDO theorem and no common squared-minor certificate. Triangle failure elsewhere is not a full-hull violation.',
            'sourceSha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (source,VERIFY)}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--check',type=Path)
    args = parser.parse_args()
    resource.setrlimit(resource.RLIMIT_CPU,(60,60))
    result = report()
    if args.check and json.loads(args.check.read_text()) != result:
        raise SystemExit('Stored certificate does not match exact recomputation')
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'interval':result['interval'],'polynomialIdentitiesVerified':True,
                      'complexPlaneCoefficientIdentitiesVerified':True,'exactParameterChecks':len(result['exactParameterChecks'])}))

if __name__ == '__main__':
    main()
