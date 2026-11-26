"""Replay fixed second-eigenvalue slices of the public four-cycle with exact arithmetic."""
from fractions import Fraction
from itertools import product
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "four-cycle-spectral"))
import spectral as s


def samples(t):
    # The other diagonal entries range over a unisolvent multi-affine grid.
    fixed = s.Q(-2) + s.Q(1, -2)*t
    for first, third, fourth in product((s.Z, s.O), repeat=3):
        yield [first, fixed, third, fourth]


def check():
    s.source_check()
    assert s.U[1][3] == s.Z and all(x != s.Z for x in s.U[1][:3])
    data = json.loads((HERE / "archive-slices.json").read_text())
    assert set(data) == {"-1", "0", "1", "-1/2", "2"}
    for parameter, entries in data.items():
        weights = {tuple(int(ch)-1 for ch in label): Fraction(value) for label, value in entries.items()}
        assert all(sorted(p) == list(range(4)) for p in weights)
        assert sum(weights.values()) == 1 and all(w > 0 for w in weights.values())
        for diagonal in samples(Fraction(parameter)):
            exact = s.e.determinant(s.matrix_sum(diagonal))
            assert exact == sum((w*s.vertex(diagonal, p) for p, w in weights.items()), s.Z)
        print(f"t={parameter}: {len(weights)} positive masses; all eight identities match")


if __name__ == "__main__":
    check()
