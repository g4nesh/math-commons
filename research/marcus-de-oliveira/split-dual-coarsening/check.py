"""Exact coverage and separating certificates for all proper spectral coarsenings."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
import importlib.util
import json

HERE = Path(__file__).resolve().parent
FINE = HERE.parent / 'split-dual-rays'
spec = importlib.util.spec_from_file_location('fine_algebra', FINE / 'algebra.py')
a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a)


def nullspace(matrix):
    red, pivots = a.rref(matrix)
    result = []
    for j in range(33):
        if j in pivots:
            continue
        v = [F(0)] * 33
        v[j] = 1
        for i, k in enumerate(pivots):
            v[k] = -red[i][j]
        result.append(a.primitive(v))
    return result


def spaces(keys, matrix):
    result = {}
    for side, indices in [('row', range(4)), ('column', range(4, 7))]:
        for i, j in combinations(indices, 2):
            grouped = {}
            for key, row in zip(keys, matrix):
                k = list(key)
                k[i] += k[j]
                del k[j]
                k = tuple(k)
                grouped.setdefault(k, [0] * 33)
                grouped[k] = [x+y for x, y in zip(grouped[k], row)]
            monomials = sorted(grouped)
            rows = [grouped[k] for k in monomials]
            _, independent = a.rref(list(map(list, zip(*rows))))
            offset = 0 if side == 'row' else 4
            name = f'{side}{i-offset}{j-offset}'
            result[name] = dict(side=side, merge=[i-offset, j-offset],
                                monomials=monomials, matrix=rows,
                                basis=independent, kernel=nullspace(rows))
    return result


def member_combination(ray, space):
    basis = [space['matrix'][j] for j in space['basis']]
    reduced, pivots = a.rref([[row[i] for row in basis] + [ray[i]] for i in range(33)])
    a.require(pivots == list(range(len(basis))), 'ray belongs to coarsened row space')
    return [[j, str(reduced[i][-1])] for i, j in enumerate(space['basis']) if reduced[i][-1]]


def validate_outcome(ray, space, outcome):
    if outcome['kind'] == 'member':
        reconstructed = [F(0)] * 33
        seen = set()
        for index, scalar in outcome['combination']:
            a.require(type(index) is int and 0 <= index < len(space['matrix']) and index not in seen,
                      'distinct valid coefficient row indices')
            seen.add(index)
            weight = F(scalar)
            a.require(weight != 0, 'nonzero rational combination entries')
            reconstructed = [x + weight*y for x, y in zip(reconstructed, space['matrix'][index])]
        a.require(reconstructed == list(ray), 'exact positive-ray row-space certificate')
    else:
        a.require(outcome['kind'] == 'outside', 'explicit outcome kind')
        index = outcome['kernel_index']
        a.require(type(index) is int and 0 <= index < len(space['kernel']), 'valid kernel witness')
        witness = space['kernel'][index]
        a.require(all(sum(x*y for x,y in zip(witness,row)) == 0 for row in space['matrix']),
                  'witness annihilates every coarsened coefficient row')
        dot = sum(x*y for x,y in zip(witness,ray))
        a.require(dot != 0 and type(outcome['dot']) is int and outcome['dot'] == dot,
                  'nonzero exact separating pairing')


def construct():
    tables, keys, matrix, kernel = a.structure()
    ss = spaces(keys, matrix)
    orbits = json.loads((FINE / 'orbits.json').read_text())['records']
    group = a.actions(tables)
    coefficient_rays = {a.primitive(row) for row in matrix}
    records = []
    for record in orbits:
        ray = record['ray']
        a.verify_ray(ray, kernel)
        outcomes = {}
        for name, space in ss.items():
            w = next((i for i,k in enumerate(space['kernel']) if sum(x*y for x,y in zip(k,ray))), None)
            if w is None:
                value = dict(kind='member', combination=member_combination(ray,space))
            else:
                value = dict(kind='outside', kernel_index=w,
                             dot=sum(x*y for x,y in zip(space['kernel'][w],ray)))
            validate_outcome(ray,space,value)
            outcomes[name] = value
        covered = any(x['kind'] == 'member' for x in outcomes.values())
        # The collection of all nine elementary mergers is physically invariant.
        images = {a.act(ray,g) for g in group}
        for image in images:
            image_covered = any(all(sum(x*y for x,y in zip(k,image)) == 0 for k in s['kernel'])
                                for s in ss.values())
            a.require(image_covered == covered, 'coverage constant on physical orbit')
        records.append(dict(id=record['id'], ray=ray, orbit_size=len(images),
                            coefficient=any(r in coefficient_rays for r in images),
                            covered=covered, outcomes=outcomes))
    return ss, records


def verify(certificate, ss, records):
    a.require(certificate['schema'] == 'mdo.split-coarsening.v1', 'certificate schema')
    a.require(len(certificate['records']) == len(records) == 62, 'all 62 physical representatives')
    for saved, actual in zip(certificate['records'],records):
        for field in ('id','ray','orbit_size','coefficient','covered'):
            a.require(saved[field] == actual[field], 'exact orbit identity and coverage')
        a.require(set(saved['outcomes']) == set(ss), 'all nine coarsenings tested')
        for name, outcome in saved['outcomes'].items():
            validate_outcome(actual['ray'],ss[name],outcome)
            a.require(outcome['kind'] == actual['outcomes'][name]['kind'], 'correct membership decision')
    return True


def run(write=False):
    ss, records = construct()
    certificate = dict(schema='mdo.split-coarsening.v1', records=records)
    if write:
        (HERE/'certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
    verify(json.loads((HERE/'certificate.json').read_text()),ss,records)
    covered = [r for r in records if r['covered']]
    new = [r for r in covered if not r['coefficient']]
    remaining = [r for r in records if not r['covered'] and not r['coefficient']]
    a.require(len(covered)==26 and len(new)==17 and len(remaining)==35, 'exact coverage totals')
    report = dict(spaces={name:dict(rank=len(s['basis']),coefficient_rows=len(s['matrix'])) for name,s in ss.items()},
                  covered_orbits=len(covered),covered_rays=sum(r['orbit_size'] for r in covered),
                  covered_noncoefficient_orbits=[r['id'] for r in new],
                  covered_noncoefficient_rays=sum(r['orbit_size'] for r in new),
                  remaining_noncoefficient_orbits=[r['id'] for r in remaining],
                  remaining_noncoefficient_rays=sum(r['orbit_size'] for r in remaining),
                  uncovered_coefficient_orbits=[r['id'] for r in records if not r['covered'] and r['coefficient']],
                  interpretation='exact linear-conic coarsening coverage; no physical negativity asserted')
    if write:
        (HERE/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    a.require(report == json.loads((HERE/'report.json').read_text()), 'saved report matches exact replay')
    return report


if __name__ == '__main__':
    import sys
    print(json.dumps(run(write='--write' in sys.argv[1:]),indent=2))
