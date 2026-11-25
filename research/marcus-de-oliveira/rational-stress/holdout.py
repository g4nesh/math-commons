#!/usr/bin/env python3
"""Separate 24-case n=7 holdout; leaves the n<=6 verifier API and report intact.

Uses its dimension-generic exact arithmetic primitives directly, with explicit
n=7 shape checks. Public certificates/verify.py intentionally still rejects n=7.
"""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import platform
import signal
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import screen as S
V = S.V


def specs():
    for repetition in range(4):
        for family in ('generic','thin-complex','clustered'):
            for mixing in ('generic','near-vertex'):
                case_id = 'n7-%s-%s-%02d' % (family,mixing,repetition)
                seed = int(hashlib.sha256(('holdout20261125:'+case_id).encode()).hexdigest()[:16],16)
                yield dict(id=case_id,n=7,family=family,mixing=mixing,repetition=repetition,seed=seed)


def analyze(witness):
    if len(witness['a']) != 7 or len(witness['b']) != 7 or len(witness['U']) != 7 or any(len(row) != 7 for row in witness['U']):
        raise ValueError('This separate holdout harness requires exactly n=7.')
    a,b = [list(map(V.parse_complex,witness[key])) for key in ('a','b')]
    U = [list(map(V.parse_complex,row)) for row in witness['U']]
    star = V.adjoint(U)
    if V.multiply(U,star) != V.identity(7) or V.multiply(star,U) != V.identity(7):
        raise V.InvalidWitness('Holdout U is not exactly unitary')
    A, B = V.diagonal(a), V.multiply(V.multiply(U,V.diagonal(b)),star)
    for matrix in (A,B):
        if V.multiply(matrix,V.adjoint(matrix)) != V.multiply(V.adjoint(matrix),matrix):
            raise ArithmeticError('Exact normality check failed')
    value = V.determinant([[A[i][j]+B[i][j] for j in range(7)] for i in range(7)])
    points = V.permutation_points(a,b)
    representatives = {}
    for permutation, point in points: representatives.setdefault(point.pair(),permutation)
    certificate = V.certify_hull(list(representatives),value.pair())
    commutator_norm = sum(S.squared_abs((a[i]-a[j])*B[i][j]) for i in range(7) for j in range(7))
    hull = certificate['hull']
    margin, edge = None,None
    if len(hull) >= 3:
        scale = max(max(p[k] for p in hull)-min(p[k] for p in hull) for k in (0,1))
        candidates = []
        for index, first in enumerate(hull):
            second = hull[(index+1)%len(hull)]
            candidates.append((V.cross(first,second,value.pair())/(sum(abs(x-y) for x,y in zip(first,second))*scale),index))
        margin, edge = min(candidates)
    membership = certificate['membership']
    record = dict(dimension=7,determinant=value.encode(),singularSum=value == V.ZERO,
                  permutationCount=len(points),uniquePointCount=len(representatives),hullDimension=min(2,len(hull)-1),membership=membership,
                  witnessSha256=S.digest(witness),noncommuting=commutator_norm != 0,commutatorNormSquared=str(commutator_norm),
                  complexUnitary=any(value.imag != 0 for row in U for value in row),
                  relativeLinfSupportMargin=str(margin) if margin is not None else None,closestEdge=edge,
                  hull=[[str(x),str(y)] for x,y in hull],exoticPsi=None,
                  classification='outside' if membership == 'outside' else 'degenerate-inside' if len(hull)<3 else 'boundary' if margin==0 else 'interior')
    if membership == 'inside':
        record['convexCombination'] = [dict(permutation=list(representatives[point]),point=list(map(str,point)),weight=str(weight)) for point,weight in certificate['weights']]
    else:
        record['separator'] = dict(normal=list(map(str,certificate['normal'])),bound=str(certificate['bound']),targetValue=str(certificate['targetValue']),gap=str(certificate['gap']))
    return record


def main():
    start_cpu,start_wall = time.process_time(),time.monotonic()
    def stop(received,_frame):
        raise S.BudgetStop(('CPU' if received == signal.SIGPROF else 'Wall')+' holdout cap reached; incomplete case omitted')
    signal.signal(signal.SIGPROF,stop); signal.signal(signal.SIGALRM,stop)
    signal.setitimer(signal.ITIMER_PROF,120); signal.setitimer(signal.ITIMER_REAL,150)
    records,reason = [],'all-24-cases-complete'
    folder = HERE/'holdout'
    (folder/'witnesses').mkdir(parents=True,exist_ok=True)
    try:
        for spec in specs():
            witness = S.make_witness(spec)
            row = analyze(witness); row['spec'] = spec
            records.append(row)
            S.write_json(folder/'witnesses'/(spec['id']+'.json'),witness)
            print(json.dumps({'case':spec['id'],'classification':row['classification'],'completed':len(records)}),flush=True)
    except S.BudgetStop as error:
        reason = str(error)
    finally:
        signal.setitimer(signal.ITIMER_PROF,0); signal.setitimer(signal.ITIMER_REAL,0)
    report = dict(schemaVersion=1,seed='holdout20261125',
                  arithmetic='exact Gaussian rationals; no numerical tolerance',
                  scope='Separate n=7 finite holdout. Does not extend the public n<=6 verifier API and proves no universal statement.',
                  budget=dict(maxCases=24,cpuSeconds=120,wallSeconds=150),
                  execution=dict(stopReason=reason,cpuSeconds=time.process_time()-start_cpu,wallSeconds=time.monotonic()-start_wall,python=platform.python_version(),platform=platform.platform()),
                  sources=dict(holdoutSha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),screenSha256=hashlib.sha256((HERE/'screen.py').read_bytes()).hexdigest(),verifierSha256=hashlib.sha256((HERE.parent/'certificates'/'verify.py').read_bytes()).hexdigest()),
                  summary=S.summarize(records),selectedWitnesses=[row['spec']['id'] for row in records],records=records)
    S.write_json(folder/'report.json',report)
    print(json.dumps(report['summary'],sort_keys=True))


if __name__ == '__main__': main()
