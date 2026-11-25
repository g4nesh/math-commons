#!/usr/bin/env python3
"""Rebuild inputs; independently check determinants and exact certificates.

Arithmetic below uses Fraction pairs, not the screen's QComplex class. The
determinant is division-free subset dynamic programming, not elimination.
This is author-side differential verification, not independent-owner review.
"""
from fractions import Fraction as F
from itertools import permutations
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import screen as S

ZERO, ONE = (F(0), F(0)), (F(1), F(0))
def add(x,y): return x[0]+y[0], x[1]+y[1]
def sub(x,y): return x[0]-y[0], x[1]-y[1]
def mul(x,y): return x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0]
def conj(x): return x[0], -x[1]
def parse(value): return F(value[0]), F(value[1])
def total(values):
    result = ZERO
    for value in values: result = add(result,value)
    return result
def cross(a,b,c): return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
def dot(a,b): return a[0]*b[0]+a[1]*b[1]
def require(condition, text):
    if not condition: raise ArithmeticError(text)


def determinant(matrix):
    n = len(matrix)
    states = {0: ONE}
    for row in range(n):
        next_states = {}
        for mask, value in states.items():
            for column in range(n):
                if mask & (1 << column): continue
                term = mul(value,matrix[row][column])
                if sum(bool(mask & (1 << k)) for k in range(column+1,n)) % 2:
                    term = (-term[0],-term[1])
                key = mask | (1 << column)
                next_states[key] = add(next_states.get(key,ZERO),term)
        states = next_states
    return states[(1 << n)-1]


def audit(witness, record):
    a, b = list(map(parse,witness['a'])), list(map(parse,witness['b']))
    U = [list(map(parse,row)) for row in witness['U']]
    n = len(a)
    require(n == len(b) == len(U) and all(len(row) == n for row in U), 'Shape mismatch')
    for i in range(n):
        for j in range(n):
            expected = ONE if i == j else ZERO
            require(total(mul(U[i][k],conj(U[j][k])) for k in range(n)) == expected, 'UU* != I')
            require(total(mul(conj(U[k][i]),U[k][j]) for k in range(n)) == expected, 'U*U != I')
    B = [[total(mul(mul(U[i][k],b[k]),conj(U[j][k])) for k in range(n)) for j in range(n)] for i in range(n)]
    matrix = [[add(B[i][j],a[i] if i == j else ZERO) for j in range(n)] for i in range(n)]
    target = determinant(matrix)
    require(target == parse(record['determinant']), 'Determinant mismatch')
    require((target == ZERO) == record['singularSum'], 'Singular flag mismatch')
    commutator = [mul(sub(a[i],a[j]),B[i][j]) for i in range(n) for j in range(n)]
    norm_squared = sum(dot(value,value) for value in commutator)
    require(norm_squared == F(record['commutatorNormSquared']), 'Commutator norm mismatch')
    require((norm_squared != 0) == record['noncommuting'], 'Noncommuting flag mismatch')
    require(any(value[1] != 0 for row in U for value in row) == record['complexUnitary'], 'Complex flag mismatch')
    if n == 4:
        def minor(rows,cols):
            z = sub(mul(U[rows[0]][cols[0]],U[rows[1]][cols[1]]),mul(U[rows[0]][cols[1]],U[rows[1]][cols[0]]))
            return dot(z,z)
        psi = (sum(dot(U[i][i],U[i][i]) for i in range(4))-minor((0,1),(0,2))-minor((0,1),(1,3))+minor((0,2),(0,1))+minor((0,2),(2,3)))/2
        require(psi == F(record['exoticPsi']),'Exotic functional mismatch')
    points = {}
    for permutation in permutations(range(n)):
        point = ONE
        for i,j in enumerate(permutation): point = mul(point,add(a[i],b[j]))
        points[tuple(j+1 for j in permutation)] = point
    unique = set(points.values())
    require(len(points) == record['permutationCount'] and len(unique) == record['uniquePointCount'], 'Point count mismatch')
    if record['membership'] == 'inside':
        weights = record['convexCombination']
        require(1 <= len(weights) <= 3, 'Invalid support size')
        require(all(F(item['weight']) >= 0 for item in weights), 'Negative weight')
        require(sum(F(item['weight']) for item in weights) == 1, 'Weight sum mismatch')
        reconstructed = ZERO
        for item in weights:
            point = points[tuple(item['permutation'])]
            require(point == parse(item['point']), 'Support point mismatch')
            reconstructed = add(reconstructed,mul(point,(F(item['weight']),F(0))))
        require(reconstructed == target, 'Convex combination mismatch')
    elif record['membership'] == 'outside':
        separator = record['separator']
        normal, bound, gap = parse(separator['normal']), F(separator['bound']), F(separator['gap'])
        require(normal != ZERO and all(dot(normal,p) <= bound for p in unique), 'Invalid separation bound')
        require(dot(normal,target)-bound == gap and gap > 0, 'Invalid positive separation')
    else: raise ArithmeticError('Unknown membership')
    hull = list(map(parse,record['hull']))
    require(hull and len(set(hull)) == len(hull) and all(p in unique for p in hull), 'Invalid hull vertices')
    dimension = 0 if len(unique) == 1 else 1 if all(cross(hull[0],hull[1],p) == 0 for p in unique) else 2
    require(dimension == record['hullDimension'], 'Hull dimension mismatch')
    if dimension == 0:
        require(len(hull) == 1, 'Invalid point hull')
    elif dimension == 1:
        require(len(hull) == 2, 'Invalid segment hull')
        for p in unique:
            require(all(min(hull[0][k],hull[1][k]) <= p[k] <= max(hull[0][k],hull[1][k]) for k in (0,1)), 'Segment omits point')
    else:
        scale = max(max(p[k] for p in unique)-min(p[k] for p in unique) for k in (0,1))
        slacks = []
        for index, first in enumerate(hull):
            second = hull[(index+1)%len(hull)]
            require(cross(first,second,hull[(index+2)%len(hull)]) > 0, 'Hull not strictly counterclockwise')
            require(all(cross(first,second,p) >= 0 for p in unique), 'Hull excludes a permutation point')
            slacks.append(cross(first,second,target)/(sum(abs(x-y) for x,y in zip(first,second))*scale))
        require(min(slacks) == F(record['relativeLinfSupportMargin']), 'Support margin mismatch')
        require(slacks[record['closestEdge']] == min(slacks), 'Closest edge mismatch')
    expected = 'outside' if record['membership'] == 'outside' else 'degenerate-inside' if dimension < 2 else 'boundary' if F(record['relativeLinfSupportMargin']) == 0 else 'interior'
    require(record['classification'] == expected, 'Classification mismatch')
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report', nargs='?', default=str(HERE/'report.json'))
    parser.add_argument('--output',default=str(HERE/'replay.json'))
    args = parser.parse_args()
    start = time.process_time()
    report_path = Path(args.report)
    report = json.loads(report_path.read_text())
    for index,row in enumerate(report['records']):
        witness = S.make_witness(row['spec'])
        require(S.digest(witness) == row['witnessSha256'],'Regenerated witness hash mismatch')
        audit(witness,row)
        if row['spec']['id'] in report['selectedWitnesses']:
            saved = json.loads((report_path.parent/'witnesses'/(row['spec']['id']+'.json')).read_text())
            require(saved == witness,'Saved witness mismatch')
        if (index+1) % 100 == 0: print('Audited %d exact cases' % (index+1),flush=True)
    require(S.summarize(report['records']) == report['summary'],'Summary mismatch')
    result = {'schemaVersion':1,'cases':len(report['records']),'mismatches':0,
              'checks':['regenerated input hash','exact unitarity','division-free determinant','commutator norm','n4 exotic functional','every permutation point','direct convex weights or separation','hull facets and margin','selected witness files','summary'],
              'cpuSeconds':time.process_time()-start,'python':platform.python_version(),
              'reportSha256':hashlib.sha256(report_path.read_bytes()).hexdigest(),
              'scope':'Author-side differential exact replay, not independent-owner review or a universal theorem.'}
    S.write_json(args.output,result)
    print(json.dumps(result,sort_keys=True))


if __name__ == '__main__': main()
