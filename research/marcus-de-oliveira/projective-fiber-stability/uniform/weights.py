"""Find and check exact shared weights for a Gaussian-rational middle and fourth eigenvalue."""
import argparse
from fractions import Fraction as F
import importlib.util
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("collision_system", HERE.parent / "check.py")
collision = importlib.util.module_from_spec(spec)
spec.loader.exec_module(collision)
sp = collision.sp


def gap_matrix(z, delta):
    rows = [("real", [{(0, 0): sp.O} for _ in collision.PERMS])]
    for subset in ((1,), (0,), (2,), (0, 2)):
        polys = [collision.translated(col, subset) if subset == (1,) else collision.pencil(col, subset, sp.Z)
                 for col in collision.initial_gaps()]
        rows.extend((part, polys) for part in ("real", "imag"))
    return collision.matrix(rows, delta, z)


def positive_law(matrix):
    """Exact Bland pivots minimize artificial mass, then maximize a common floor."""
    count = len(matrix[0])
    height = len(matrix)
    # m_j = excess_j + floor. Artificial identity columns give a feasible initial basis.
    tableau = [row + [sum(row)] + [F(i == k) for k in range(height)] + [F(i == 0)]
               for i, row in enumerate(matrix)]
    costs = [(0, 0)] * count + [(0, -1)] + [(1, 0)] * height
    basis = list(range(count + 1, count + 1 + height))
    seen = set()
    while tuple(basis) not in seen:
        seen.add(tuple(basis))
        entering = None
        for j, cost in enumerate(costs):
            if j in basis:
                continue
            reduced = tuple(cost[h] - sum(costs[b][h] * tableau[i][j] for i, b in enumerate(basis))
                            for h in (0, 1))
            if reduced < (0, 0):
                entering = j
                break
        if entering is None:
            break
        candidates = [(row[-1] / row[entering], basis[i], i) for i, row in enumerate(tableau)
                      if row[entering] > 0]
        if not candidates:
            raise ArithmeticError("Unbounded normalized probability tableau")
        _, _, leaving = min(candidates)
        pivot = [v / tableau[leaving][entering] for v in tableau[leaving]]
        for i, row in enumerate(tableau):
            if i != leaving and row[entering]:
                tableau[i] = [a - row[entering] * b for a, b in zip(row, pivot)]
        tableau[leaving] = pivot
        basis[leaving] = entering
    else:
        raise ArithmeticError("Repeated exact simplex basis")

    solution = [F(0)] * len(costs)
    for i, j in enumerate(basis):
        solution[j] = tableau[i][-1]
    if any(solution[count + 1:]) or solution[count] <= 0:
        raise ValueError("No strictly positive solution to these regularized coefficient equations")
    masses = [solution[count] + v for v in solution[:count]]
    if min(masses) <= 0 or any(sum(a*m for a, m in zip(row, masses)) != int(i == 0)
                               for i, row in enumerate(matrix)):
        raise ArithmeticError("Constructed law failed direct exact substitution")
    return masses


def weights_at(z, delta):
    return positive_law(gap_matrix(z, delta))


def replay():
    # The same probabilities must work at unrelated x,y and at repeated eigenvalues.
    for z, delta in [(sp.Q(F(3, 5), F(7, 10)), sp.Q(F(1, 2048), -F(1, 4096))),
                     (sp.I, sp.Q(F(1, 1024))), (sp.Z, sp.Z), (sp.Q(120, -37), sp.Q(-F(1, 8192)))]:
        masses = weights_at(z, delta)
        for x, y in ((sp.Q(-2, 3), sp.Q(4, -1)), (z, z)):
            a = [x, z, y, delta]
            actual = sp.e.determinant(sp.matrix_sum(a))
            weighted = sum((m * sp.vertex(a, p) for m, p in zip(masses, collision.PERMS)), sp.Z)
            if weighted != actual:
                raise ArithmeticError("Original dense determinant disagrees with the common law")
    print("Exact positive laws verified against original determinants on four fibers.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("coordinates", nargs="*", type=F, metavar="RATIONAL")
    tokens = sys.argv[1:]
    if tokens in (["-h"], ["--help"]):
        parser.print_help()
        sys.exit(0)
    if tokens[:1] == ["--"]:
        tokens = tokens[1:]
    # argparse otherwise treats negative fractions such as -1/4096 as options.
    args = parser.parse_args(["--", *tokens])
    if not args.coordinates:
        replay()
    elif len(args.coordinates) != 4:
        parser.error("Supply Re(z) Im(z) Re(delta) Im(delta), or no arguments for the replay.")
    else:
        coords = args.coordinates
        law = weights_at(sp.Q(*coords[:2]), sp.Q(*coords[2:]))
        print(json.dumps({"".join(str(i+1) for i in p): str(m) for p, m in zip(collision.PERMS, law)}, indent=2))
