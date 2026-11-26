"""Replay the two bounded exact coefficient discovery scans; no positivity inference."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import json
import platform
import random
import time
import verify

HERE = Path(__file__).resolve().parent
r = verify.r
SEED = 20261125


def run(cases, zero_u11):
    rng = random.Random(SEED)
    start = time.process_time()

    def pair():
        n = rng.randint(2, 30)
        result = (F(n*n-1, n*n+1), F(2*n, n*n+1))
        return result if rng.randrange(2) else result[::-1]

    def phase():
        a, b = pair()
        return r.z(a if rng.randrange(2) else -a, b if rng.randrange(2) else -b)

    minima = [None]*4
    negatives = [0]*4
    worst = None
    for case in range(cases):
        p, q = pair()
        co, si = (F(0), F(1)) if zero_u11 else pair()
        ce, te = pair()
        ra, sa = pair()
        ar, rr = pair()
        U = r.boundary(p, q, co, si, ce, te, r.scale(ra, phase()), r.scale(sa, phase()),
                       r.z(F(3, 5)), r.z(F(4, 5)), r.scale(ar, phase()), rr)
        data = verify.scalars(U)  # Includes exact unitarity validation.
        coefficients = data['coefficients']
        values = [2*data['k']+data['h'], data['k'], coefficients[2], coefficients[3]]
        for j, value in enumerate(values):
            if minima[j] is None or value < minima[j]:
                minima[j] = value
            negatives[j] += value < 0
        if coefficients[2] < 0 and (worst is None or coefficients[2] < F(worst['coefficients'][2])):
            worst = {'case': case, 'U': [[[str(re), str(im)] for re, im in row] for row in U],
                     'coefficients': list(map(str, coefficients)),
                     **{key: str(data[key]) for key in ('k', 'h', 'd', 'e')},
                     'schema': 'exact-gaussian-rational-matrix-v1',
                     'claim': 'Negative cubic t^2 coefficient only; not a counterexample to star.'}
    name = 'u11-negative-c2.json' if zero_u11 else 'negative-c2.json'
    saved = json.loads((HERE / name).read_text())
    assert worst == saved, 'Discovery replay differs from the saved exact strongest fixture.'
    return {'cases': cases, 'seed': SEED, 'u11_zero': zero_u11,
            'minima_labels': ['2k+h', 'k', 'c2', 'c3'],
            'minima_exact': list(map(str, minima)), 'negative_counts': negatives,
            'saved_strongest_fixture_reproduced': name, 'python': platform.python_version(),
            'cpu_seconds': time.process_time()-start,
            'scope': 'Finite parametrized exact coefficient scan; no global sign conclusion.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-reports', action='store_true')
    args = parser.parse_args()
    reports = [run(1500, False), run(1000, True)]
    if args.write_reports:
        (HERE / 'scan-report.json').write_text(json.dumps(reports, indent=2)+'\n')
    print(json.dumps(reports, indent=2))
