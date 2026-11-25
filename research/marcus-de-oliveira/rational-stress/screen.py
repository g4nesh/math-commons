#!/usr/bin/env python3
"""Bounded exact MDO instances, constructed in Q(i), never rounded from floats."""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import platform
import random
import signal
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'certificates'))
import verify as V

SEED = 20261125
FAMILIES = ('generic', 'thin-complex', 'clustered', 'singular', 'two-value', 'hermitian')
MIXING = ('generic', 'near-vertex', 'exotic-centered')
Q = V.QComplex


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=True).encode()


def digest(value):
    return hashlib.sha256(canonical(value)).hexdigest()


def write_json(path, value):
    with Path(path).open('w', encoding='utf8', newline='\n') as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write('\n')


def case_specs(seed=SEED):
    """Round-robin dimensions so a budget stop does not omit n=6."""
    for repetition in range(20):
        for family in FAMILIES:
            for mixing in MIXING:
                for n in (4, 5, 6):
                    if n == 6 and repetition >= 10:
                        continue
                    case_id = 'n%d-%s-%s-%02d' % (n, family, mixing, repetition)
                    case_seed = int(hashlib.sha256(('%d:%s' % (seed, case_id)).encode()).hexdigest()[:16], 16)
                    yield dict(id=case_id, n=n, family=family, mixing=mixing, repetition=repetition, seed=case_seed)


def spectra(spec, rng):
    n, family = spec['n'], spec['family']
    randq = lambda: Q(rng.randrange(-5, 6), rng.randrange(-5, 6))
    a = [Q(2*i-n, rng.randrange(-4, 5)) for i in range(n)]
    b = [Q(3*i-n+1, rng.randrange(-5, 6)) for i in range(n)]
    if family == 'thin-complex':
        epsilon = F(1, 2 ** (8 + 2*(spec['repetition'] % 6)))
        a = [Q(2*i+1, epsilon*rng.randrange(-4, 5)) for i in range(n)]
        b = [Q(3*i+2, epsilon*rng.randrange(-5, 6)) for i in range(n)]
    elif family == 'clustered':
        epsilon = F(1, 2 ** (4 + 2*(spec['repetition'] % 6)))
        center_a, center_b = randq(), randq()
        a = [center_a + epsilon*value for value in a]
        b = [center_b + epsilon*value for value in b]
    elif family == 'singular':
        # U fixes coordinate zero; the remaining block is still noncommuting.
        b[0] = -a[0]
    elif family == 'two-value':
        left, right = Q(-2, 1), Q(3, -2)
        b = [left if i < n//2 else right for i in range(n)]
    elif family == 'hermitian':
        a = [Q(2*i-n) for i in range(n)]
        b = [Q(3*i-n+1) for i in range(n)]
    if spec['mixing'] == 'near-vertex' and family != 'singular':
        # Lexicographically minimal sigma point is an exposed hull vertex.
        permutation, _ = min(V.permutation_points(a, b), key=lambda item: item[1].pair())
        b = [b[index-1] for index in permutation]
    return a, b


def rotate_rows(unitary, i, j, t, phase):
    """Left multiply by [[c,-s*conj(z)],[s*z,c]], |z|=1, t in Q."""
    c, s = (1-t*t)/(1+t*t), 2*t/(1+t*t)
    first, second = unitary[i][:], unitary[j][:]
    unitary[i] = [c*x-s*phase.conjugate()*y for x, y in zip(first, second)]
    unitary[j] = [s*phase*x+c*y for x, y in zip(first, second)]


def make_witness(spec):
    rng = random.Random(spec['seed'])
    n = spec['n']
    a, b = spectra(spec, rng)
    unitary = V.identity(n)
    start = 1 if spec['family'] == 'singular' else 0
    if spec['mixing'] == 'exotic-centered' and n-start >= 4:
        rows = ((1,4,4,4), (-4,1,4,-4), (-4,-4,1,4), (-4,4,-4,1))
        for i in range(4):
            for j in range(4):
                unitary[start+i][start+j] = Q(F(rows[i][j], 7))
    indices = list(range(start, n))
    rng.shuffle(indices)
    pairs = list(zip(indices, indices[1:]))
    pairs += [tuple(rng.sample(indices, 2)) for _ in range(n//2)]
    phases = (Q(0,1), Q(F(3,5),F(4,5)), Q(F(5,13),F(-12,13)))
    for index, (i, j) in enumerate(pairs):
        if spec['mixing'] == 'generic':
            t = rng.choice((F(1,2), F(2,3), F(1,3), F(3,4)))
        else:
            t = F(1, 2 ** ((4,8,10)[spec['repetition'] % 3]))
        if rng.randrange(2):
            t = -t
        rotate_rows(unitary, i, j, t, phases[(index+rng.randrange(3)) % 3])
    return {'schemaVersion': 1, 'name': spec['id'],
            'notes': 'Exact rational construction; finite instance, not a universal proof.',
            'a': [value.encode() for value in a], 'b': [value.encode() for value in b],
            'U': [[value.encode() for value in row] for row in unitary]}


def squared_abs(value):
    return value.real*value.real + value.imag*value.imag


def exotic_psi(unitary):
    """Zheng functional transcribed from minors/moc/crosscheck_zheng.py, n=4 only."""
    minor = lambda rows, cols: squared_abs(unitary[rows[0]][cols[0]]*unitary[rows[1]][cols[1]] - unitary[rows[0]][cols[1]]*unitary[rows[1]][cols[0]])
    return (sum(squared_abs(unitary[i][i]) for i in range(4))
            - minor((0,1),(0,2)) - minor((0,1),(1,3))
            + minor((0,2),(0,1)) + minor((0,2),(2,3))) / 2


def analyze(witness):
    result = V.verify_witness(witness)
    a, b, unitary, _ = V.parse_witness(witness)
    matrix_b = V.multiply(V.multiply(unitary, V.diagonal(b)), V.adjoint(unitary))
    commutator = [[(a[i]-a[j])*matrix_b[i][j] for j in range(len(a))] for i in range(len(a))]
    norm_squared = sum(squared_abs(value) for row in commutator for value in row)
    hull = [tuple(map(F, point)) for point in result['hull']]
    target = tuple(map(F, result['determinant']))
    margin, supporting = None, None
    if len(hull) >= 3:
        scale = max(max(p[k] for p in hull)-min(p[k] for p in hull) for k in (0,1))
        candidates = []
        for index, first in enumerate(hull):
            second = hull[(index+1) % len(hull)]
            dx, dy = second[0]-first[0], second[1]-first[1]
            candidates.append((V.cross(first, second, target)/((abs(dx)+abs(dy))*scale), index))
        margin, supporting = min(candidates)
    classification = 'outside' if result['membership'] == 'outside' else ('degenerate-inside' if len(hull) < 3 else 'boundary' if margin == 0 else 'interior')
    record = {key: result[key] for key in ('dimension','determinant','singularSum','permutationCount','uniquePointCount','hullDimension','membership')}
    record.update(witnessSha256=digest(witness), classification=classification,
                  noncommuting=norm_squared != 0, commutatorNormSquared=str(norm_squared),
                  complexUnitary=any(value.imag for row in unitary for value in row),
                  relativeLinfSupportMargin=str(margin) if margin is not None else None,
                  closestEdge=supporting, hull=result['hull'],
                  exoticPsi=str(exotic_psi(unitary)) if len(a) == 4 else None)
    # bool conversion avoids Fraction leaking into JSON.
    record['complexUnitary'] = bool(record['complexUnitary'])
    if result['membership'] == 'inside':
        record['convexCombination'] = result['convexCombination']
    else:
        record['separator'] = result['separator']
    return record


def summarize(records):
    groups = {}
    for key in ('n', 'family', 'mixing'):
        groups[key] = {str(value): dict(Counter(row['classification'] for row in records if row['spec'][key] == value))
                       for value in sorted({row['spec'][key] for row in records})}
    return {'cases': len(records), 'noncommuting': sum(row['noncommuting'] for row in records),
            'complexUnitary': sum(row['complexUnitary'] for row in records), 'singularSums': sum(row['singularSum'] for row in records),
            'classification': dict(Counter(row['classification'] for row in records)), 'groups': groups,
            'n4NegativeExoticPsi': sum(row['exoticPsi'] is not None and F(row['exoticPsi']) < 0 for row in records),
            'generalComplexCases': sum(row['spec']['family'] in ('generic','thin-complex','clustered') for row in records),
            'structuredControls': sum(row['spec']['family'] in ('singular','two-value','hermitian') for row in records),
            'permutationPointsChecked': sum(row['permutationCount'] for row in records)}


class BudgetStop(Exception):
    pass


def run(args):
    started_cpu, started_wall = time.process_time(), time.monotonic()
    def stop(received, _frame):
        raise BudgetStop(('CPU' if received == signal.SIGPROF else 'Wall-time') + ' cap reached during a case; incomplete case omitted.')
    signal.signal(signal.SIGALRM, stop)
    signal.signal(signal.SIGPROF, stop)
    signal.setitimer(signal.ITIMER_REAL, args.wall_seconds)
    signal.setitimer(signal.ITIMER_PROF, args.cpu_seconds)
    records, reason = [], 'case-limit-complete'
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    selected = output / 'witnesses'
    selected.mkdir(exist_ok=True)
    try:
        for spec in case_specs(args.seed):
            if len(records) >= args.max_cases:
                break
            if time.process_time()-started_cpu >= args.cpu_seconds:
                reason = 'cpu-budget-reached'; break
            witness = make_witness(spec)
            row = analyze(witness)
            row['spec'] = spec
            records.append(row)
            if len(records) % 30 == 0:
                print(json.dumps({'cases':len(records),'cpuSeconds':round(time.process_time()-started_cpu,3),'last':spec['id']}), flush=True)
    except BudgetStop as error:
        reason = str(error)
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.setitimer(signal.ITIMER_PROF, 0)
    # Full input files for every boundary/outside case and the two tightest
    # positive-margin cases per dimension/mixing. All other inputs regenerate.
    chosen = {row['spec']['id'] for row in records if row['classification'] in ('outside','boundary')}
    for n in (4,5,6):
        for mixing in MIXING:
            eligible = [row for row in records if row['spec']['n'] == n and row['spec']['mixing'] == mixing and row['relativeLinfSupportMargin'] is not None and F(row['relativeLinfSupportMargin']) > 0]
            chosen.update(row['spec']['id'] for row in sorted(eligible, key=lambda row:F(row['relativeLinfSupportMargin']))[:2])
    for row in records:
        if row['spec']['id'] in chosen:
            write_json(selected / (row['spec']['id']+'.json'), make_witness(row['spec']))
    report = {'schemaVersion':1, 'seed':args.seed, 'arithmetic':'exact Gaussian rationals; no numerical tolerance',
              'scope':'Finite constructed instances only; no universal proof or counterexample claimed without independent review.',
              'budget':{'maxCases':args.max_cases,'cpuSeconds':args.cpu_seconds,'wallSeconds':args.wall_seconds},
              'execution':{'stopReason':reason,'cpuSeconds':time.process_time()-started_cpu,'wallSeconds':time.monotonic()-started_wall,
                           'python':platform.python_version(),'platform':platform.platform()},
              'sources':{'screenSha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                         'verifierSha256':hashlib.sha256((HERE.parent/'certificates'/'verify.py').read_bytes()).hexdigest()},
              'summary':summarize(records), 'selectedWitnesses':sorted(chosen), 'records':records}
    write_json(output/'report.json', report)
    print(json.dumps(report['summary'],sort_keys=True), flush=True)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed',type=int,default=SEED)
    parser.add_argument('--max-cases',type=int,default=900)
    parser.add_argument('--cpu-seconds',type=int,default=600)
    parser.add_argument('--wall-seconds',type=int,default=720)
    parser.add_argument('--output',default=str(HERE))
    args = parser.parse_args()
    if not 1 <= args.max_cases <= 900 or not 1 <= args.cpu_seconds <= 600 or not 1 <= args.wall_seconds <= 720:
        parser.error('Bounds:1..900 cases,1..600 CPU seconds,1..720 wall seconds.')
    run(args)


if __name__ == '__main__':
    main()
