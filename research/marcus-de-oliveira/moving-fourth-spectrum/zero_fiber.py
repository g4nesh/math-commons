"""Exact obstruction to weights uniform over three eigenvalues near a₃ = 0."""

import json
import sys
from fractions import Fraction
from itertools import permutations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'four-cycle-spectral'))
from fiber_certificate import target
from spectral import Q, Z, source_check, vertex


def replay(certificate):
    source_check()
    pairings = list(permutations(range(4)))
    samples = (Q(0), Q(1), Q(0, 1))
    gaps = {p: [Fraction(0) for _ in samples] for p in pairings}
    seen = set()
    for probe in certificate['probes']:
        bits = tuple(probe['triple'])
        if len(bits) != 3 or bits in seen or any(type(x) is not int or x not in (0, 1) for x in bits):
            raise ValueError("Expected distinct binary triples")
        seen.add(bits)
        masses = {}
        for term in probe['pointwiseWeights']:
            p, w = tuple(term['p']), Fraction(term['w'])
            if p not in gaps or p in masses or w <= 0:
                raise ValueError("Invalid pointwise permutation or weight")
            masses[p] = w
        if sum(masses.values()) != 1:
            raise ValueError("Pointwise weights must sum to one")
        x, y = Fraction(probe['real']), Fraction(probe['imag'])
        for i, z in enumerate(samples):
            diagonal = [Q(bit) for bit in bits] + [z]
            actual = target(diagonal)
            vertices = {p: vertex(diagonal, p) for p in pairings}
            # These masses work for the individual triple and every complex z,
            # since both sides are holomorphic affine in z and agree at 0 and 1.
            if sum((w * vertices[p] for p, w in masses.items()), Z) != actual:
                raise ValueError(f"Individual hull certificate failed at {bits}, {z}")
            for p, point in vertices.items():
                delta = actual - point
                gaps[p][i] += x * delta.real + y * delta.imag

    # Each weighted real part of an affine complex function is affine in x,y.
    # Its minimum on [-r,r]² is c-r(|dx|+|dy|). Take all 24 inequalities.
    if not seen or any(values[0] <= 0 for values in gaps.values()):
        raise ValueError("The fixed functional must separate every permutation at zero")
    radii = []
    for constant, at_one, at_i in gaps.values():
        variation = abs(at_one - constant) + abs(at_i - constant)
        if variation:
            radii.append(constant / variation)
    if not radii:
        raise ValueError("Expected a finite boundary for this functional")
    return {
        'triples': len(seen),
        'permutations': len(pairings),
        'minimumGapAtZero': str(min(v[0] for v in gaps.values())),
        'openSquareRadius': str(min(radii)),
    }


if __name__ == '__main__':
    print(json.dumps(replay(json.loads(Path(__file__).with_suffix('.json').read_text())), indent=2))
