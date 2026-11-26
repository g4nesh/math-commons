#!/usr/bin/env python3
"""Solver-free exact classification of the archived fourth-eigenvalue segment."""
from fractions import Fraction as F
from itertools import permutations, product
from pathlib import Path
import importlib.util
import json
import sys

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('boundary_moving', HERE.parent / 'moving-fourth-spectrum' / 'moving.py')
m = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = m
spec.loader.exec_module(m)
P = tuple(permutations(range(4)))
ALPHA = tuple(map(F, (-1400, 8398, 9419, 1808)))
BETA = tuple(map(F, (-36, 183, 212, 41)))
ALPHA_BRACKET = (F(143, 1000), F(144, 1000))
BETA_BRACKET = (F(164, 1000), F(165, 1000))
SWITCH = F(97, 100)


def system(t):
    return m.coefficient_system(m.s.Q(t, -t), P)


def require(ok, message):
    if not ok:
        raise ValueError(message)


def polynomial(p):
    require(isinstance(p, list) and bool(p) and all(type(v) is int for v in p), 'Expected integer polynomial')
    require(len(p) <= 13, 'Polynomial exceeds Cramer degree bound')
    return tuple(map(F, p))


def positive(p, bounds):
    controls = m.bernstein(p, *bounds)
    require(min(controls) > 0, 'Bernstein strict positivity failed')
    return min(controls)


def support_indices(support):
    require(isinstance(support, list) and len(support) == 12, 'Expected twelve support permutations')
    require(all(isinstance(p, list) and len(p) == 4 and all(type(v) is int for v in p) and sorted(p) == list(range(4)) for p in support), 'Invalid permutation')
    indices = [P.index(tuple(p)) for p in support]
    require(len(set(indices)) == 12, 'Duplicate support permutation')
    return indices


def cleared_identities(base, A0, b0, A1, b1):
    ids = support_indices(base['support'])
    den = polynomial(base['denominator'])
    require(len(base['numerators']) == 12, 'Wrong numerator count')
    nums = [polynomial(p) for p in base['numerators']]
    for i in range(17):
        left = (F(0),)
        for j, p in zip(ids, nums):
            left = m.padd(left, m.pmul((A0[i][j], A1[i][j] - A0[i][j]), p))
        require(left == m.pmul((b0[i], b1[i] - b0[i]), den), 'Cleared coefficient identity failed')
    return den, nums


def verify_positive(data):
    require(tuple(data['alphaPolynomial']) == ALPHA and tuple(data['betaPolynomial']) == BETA, 'Root polynomials changed')
    require(data['alphaIsolation'] == list(map(str, ALPHA_BRACKET)) and data['betaIsolation'] == list(map(str, BETA_BRACKET)), 'Root brackets changed')
    require(data['switch'] == str(SWITCH), 'Switch changed')
    # Every nonconstant coefficient is positive: both cubics are strictly
    # increasing on [0,1]. Opposite endpoint signs isolate unique roots.
    for p, bracket in ((ALPHA, ALPHA_BRACKET), (BETA, BETA_BRACKET)):
        require(all(v > 0 for v in p[1:]), 'Root monotonicity failed')
        require(m.evaluate(p, bracket[0]) < 0 < m.evaluate(p, bracket[1]), 'Root isolation failed')
    require(ALPHA_BRACKET[1] < BETA_BRACKET[0] < SWITCH, 'Intervals are unordered')
    A0, b0 = system(F(0)); A1, b1 = system(F(1))
    require(len(data['bases']) == 3, 'Expected three bases')
    checked = [cleared_identities(base, A0, b0, A1, b1) for base in data['bases']]
    bounds = ((ALPHA_BRACKET[0], BETA_BRACKET[1]), (BETA_BRACKET[0], SWITCH), (SWITCH, F(1)))
    minimum_controls = []
    for k, ((den, nums), interval) in enumerate(zip(checked, bounds)):
        values = [positive(den, interval)]
        for j, p in enumerate(nums):
            if k == 0 and j == 0:
                require(p == tuple(-729 * v for v in BETA), 'First transition factor failed')
            elif k == 0 and j == 1:
                require(p == tuple(16 * v for v in ALPHA), 'Boundary factor failed')
            elif k == 1 and j == 0:
                quotient, remainder = m.pdiv(p, BETA)
                require(remainder == (F(0),), 'Second transition factor failed')
                values.append(positive(quotient, interval))
            else:
                values.append(positive(p, interval))
        minimum_controls.append(str(min(values)))
    # Direct determinant checks are additional regression evidence. The proof
    # for arbitrary complex triples is the polynomial coefficient identity.
    for t in (F(3,20), F(1,4), SWITCH, F(1)):
        weights = weights_at(t, data)
        for triple in product((0,1), repeat=3):
            a = [m.s.Q(bit) for bit in triple] + [m.s.Q(t,-t)]
            target = m.s.leibniz(m.s.matrix_sum(a))
            require(target == m.s.e.determinant(m.s.matrix_sum(a)), 'Independent determinant disagreement')
            require(target == sum((w*m.s.vertex(a,p) for p,w in weights), m.s.Z), 'Direct corner reconstruction failed')
    return {'closedFeasibleSet':'[alpha,1]', 'alphaPolynomial':list(map(int,ALPHA)), 'alphaIsolation':list(map(str,ALPHA_BRACKET)),
            'supportIntervals':['[alpha,beta]','[beta,97/100]','[97/100,1]'], 'realPolynomialIdentities':51,
            'strictBernsteinControlMinima':minimum_controls,'betaBothFormulasFeasible':True,'directCornerChecks':32,
            'algebraicEndpointMasses':{'alphaPositive':11,'alphaZero':1,'betaPositive':11,'betaZero':1}}


def weights_at(t, data=None):
    if isinstance(t,bool) or not isinstance(t,(int,F)):
        raise TypeError('Use an exact integer or Fraction parameter')
    t = F(t)
    if not 0 <= t <= 1 or m.evaluate(ALPHA,t) < 0:
        raise ValueError('Parameter is outside the certified feasible interval')
    data = data if data is not None else json.loads((HERE/'positive.json').read_text())
    index = 0 if m.evaluate(BETA,t) <= 0 else 1 if t <= SWITCH else 2
    base = data['bases'][index]
    den = m.evaluate(base['denominator'],t)
    require(den > 0, 'Nonpositive denominator')
    weights = tuple((tuple(p),m.evaluate(n,t)/den) for p,n in zip(base['support'],base['numerators']))
    require(all(w>=0 for _,w in weights) and sum(w for _,w in weights)==1, 'Invalid convex weights')
    return weights


def verify_dual(data):
    ids = support_indices(data['support']); rows = data['basisRows']; column = data['replacedColumn']
    require(len(rows)==12 and len(set(rows))==12 and all(type(i)is int and 0<=i<17 for i in rows), 'Invalid dual rows')
    require(type(column)is int and 0<=column<12, 'Invalid replaced column')
    require(len(data['columnPolynomials'])==24, 'Expected all24 permutation columns')
    polynomials = [polynomial(p)for p in data['columnPolynomials']]+[polynomial(data['targetPolynomial'])]
    require(all(isinstance(x,str)for x in data['commonFactor']), 'Common factor requires exact rational strings')
    common = tuple(F(x)for x in data['commonFactor'])
    require(common == (F(0),F(0),F(0),F(50176,5),F(12544,5)), 'Wrong common scalar factor')
    require(all(len(m.pmul(common,p))<=13 for p in polynomials), 'Normalized polynomial exceeds the determinant degree bound')
    # det(C with one column replaced) is linear in that column. Every entry
    # is affine in s, hence each determinant has degree at most12. Thirteen
    # exact values prove the recorded polynomial identities everywhere.
    for t in range(13):
        A,b = system(F(t)); C = [[A[r][j]for j in ids]for r in rows]
        columns = [[A[r][j]for r in rows]for j in range(24)]+[[b[r]for r in rows]]
        for values,p in zip(columns,polynomials):
            actual = m.real_det([row[:column]+[v]+row[column+1:]for row,v in zip(C,values)])
            require(actual == m.evaluate(common,t)*m.evaluate(p,t), 'Cofactor polynomial identity failed')
    controls = [m.bernstein(p,F(0),ALPHA_BRACKET[1])for p in polynomials[:-1]]
    require(all(min(c)>=0 for c in controls), 'Dual permutation sign failed')
    require(polynomials[-1] == m.pmul((32,96,32),ALPHA), 'Dual target factor failed')
    # common=(12544/5)*s^3*(s+4)>0 for s>0. Thus the normalized
    # cofactor functional is nonnegative on every vertex column and negative
    # on the target for 0<s<alpha. The singular endpoint uses a separate dual.
    require(len(data['zeroDual'])==17 and all(isinstance(x,str)for x in data['zeroDual']), 'Invalid zero-endpoint dual')
    y = list(map(F,data['zeroDual'])); A,b = system(F(0))
    require(all(sum(y[i]*A[i][j]for i in range(17))<=0 for j in range(24)), 'Zero-endpoint column sign failed')
    gap = sum(x*v for x,v in zip(y,b));require(gap>0,'Zero-endpoint target gap failed')
    return {'excludedSet':'[0,alpha)','allPermutationColumns':24,'cofactorPolynomialDegreeBound':12,
            'exactInterpolationValues':13,'normalizingFactor':'(12544/5)*s^3*(s+4)',
            'targetFactor':'32*(1+3*s+s^2)*(1808*s^3+9419*s^2+8398*s-1400)',
            'zeroEndpointGap':str(gap),'minimumColumnBernsteinControl':str(min(min(c)for c in controls))}


def verify():
    m.s.source_check()
    require(m.s.e.multiply(m.s.U,m.s.e.adjoint(m.s.U))==m.s.e.identity(4),'Archived U is not exactly unitary')
    require(m.s.cross(m.s.B_SPECTRUM[0]-m.s.B_SPECTRUM[1],m.s.B_SPECTRUM[3]-m.s.B_SPECTRUM[1])==-5,'Fixed spectrum is no longer the archived noncollinear b')
    positive_data = json.loads((HERE/'positive.json').read_text());dual_data=json.loads((HERE/'dual.json').read_text())
    return {'parameter':'a3=s*(1-i), 0<=s<=1; fixed U,b and arbitrary complex a0,a1,a2',
            'positive':verify_positive(positive_data),'negative':verify_dual(dual_data),
            'scope':'Exact existence classification for one distribution uniform over three free eigenvalues. Not a pointwise MDO counterexample or a general theorem.'}


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report',type=Path,help='Write the deterministic exact replay summary')
    parser.add_argument('--weights',type=F,help='Print convex weights at a feasible rational s')
    args=parser.parse_args()
    result=verify()
    if args.report:args.report.write_text(json.dumps(result,indent=2)+'\n')
    if args.weights is not None:result['weights']=[{'columnToRow':list(p),'mass':str(w)}for p,w in weights_at(args.weights)]
    print(json.dumps(result,indent=2))
