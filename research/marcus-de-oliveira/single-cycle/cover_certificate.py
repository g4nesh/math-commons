#!/usr/bin/env python3
"""Standalone exact four-cycle obstruction checker; Python standard library only.

Uses Fraction matrices directly, not the other research arithmetic modules or
any optimizer. Enumerates all 24 permutations and uses only 1x1/2x2 minors.
"""
from fractions import Fraction as F
from itertools import permutations
import argparse
import json
from pathlib import Path

COVER = (
 ((0,), (0,)), ((1,), (3,)), ((2,), (1,)), ((3,), (2,)),
 ((0,1),(0,2)), ((0,1),(0,3)), ((0,1),(1,3)),
 ((0,2),(0,1)), ((0,2),(0,3)), ((0,2),(1,2)),
 ((0,3),(0,1)), ((0,3),(0,2)), ((0,3),(2,3)),
)


def identity(n):
    return [[F(int(i == j)) for j in range(n)] for i in range(n)]


def multiply(a, b):
    return [[sum(a[i][k]*b[k][j] for k in range(len(a))) for j in range(len(a))] for i in range(len(a))]


def matrix():
    u, c, s = identity(4), F(3,5), F(4,5)
    for i,j in ((0,1),(1,2),(2,3),(0,3)):
        gate = identity(4)
        gate[i][i],gate[i][j],gate[j][i],gate[j][j] = c,s,-s,c
        u = multiply(gate,u)
    return u


def check_cover(cover):
    if len(set(cover)) != len(cover):
        raise ValueError('Repeated cover coordinate')
    rows = []
    for p in permutations(range(4)):
        hits = sum(tuple(sorted(p[j] for j in t)) == s for s,t in cover)
        if hits < 1:
            raise ValueError('The alleged cover misses a permutation')
        rows.append({'columnToRow':list(p),'hits':hits})
    return rows


def certify():
    u = matrix()
    transpose = [list(row) for row in zip(*u)]
    assert multiply(u,transpose) == identity(4) == multiply(transpose,u)
    permutation_checks = check_cover(COVER)
    terms = []
    for s,t in COVER:
        minor = u[s[0]][t[0]] if len(s) == 1 else u[s[0]][t[0]]*u[s[1]][t[1]]-u[s[0]][t[1]]*u[s[1]][t[0]]
        terms.append({'rows':list(s),'columns':list(t),'squaredMinor':str(minor*minor)})
    value = sum(F(term['squaredMinor']) for term in terms)
    assert value == F(14953,15625) and 1-value == F(672,15625)
    return {'schemaVersion':1,'claim':'No probability distribution on S4 matches every squared minor of this exact four-cycle unitary.',
            'scope':'Obstruction to common-minor certificates, not to the Marcus-de Oliveira determinant conjecture.',
            'matrixConvention':'R14 R34 R23 R12 with blocks [[3/5,4/5],[-4/5,3/5]]; zero-based supports in data',
            'U':[[str(z) for z in row] for row in u],
            'cover':terms,'permutationChecks':permutation_checks,
            'permutationLowerBound':'1','unitaryCoverValue':str(value),'strictGap':str(1-value)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--check',type=Path)
    args = parser.parse_args()
    result = certify()
    if args.check and json.loads(args.check.read_text()) != result:
        raise SystemExit('Stored certificate differs from exact independent recomputation')
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'permutationsChecked':24,'coverCoordinates':13,
                      'permutationLowerBound':result['permutationLowerBound'],
                      'unitaryCoverValue':result['unitaryCoverValue'],'strictGap':result['strictGap'],
                      'exactObstructionVerified':True}))

if __name__ == '__main__':
    main()
