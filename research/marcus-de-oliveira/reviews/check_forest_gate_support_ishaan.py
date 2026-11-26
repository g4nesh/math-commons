"""Independent checks for mdo-forest-gate-support (formal/ForestGateSupport.lean).

Exact Gaussian-integer arithmetic, standard library only. Determinants use the
Leibniz permutation sum, not the Laplace recursion of the Lean file.
  A. Single gate (identity or arbitrary diagonal outside the 2x2 block): every
     nonzero square minor det M[T,S] has T == S or T == swap(S).
  B. Sharpness: for a generic identity-elsewhere gate every stay or swap is nonzero.
  C. Forest schedules: at most one subset history with nonzero minor product
     between fixed endpoints; the sum over histories equals the compound entry
     of the product (Cauchy-Binet; a consistency check, not claimed by the artifact).
  D. Controls: a triangle and a repeated edge each give two nonzero histories.
With --lean, also compile the reviewed modules and a probe (non-vacuity examples
and #eval of minorDet compared with Leibniz determinants).
"""
import itertools, os, random, re, subprocess, sys, tempfile
from pathlib import Path

SEED, ZERO, ONE = 20261125, (0, 0), (1, 0)
rng = random.Random(SEED)
add = lambda x, y: (x[0] + y[0], x[1] + y[1])
mul = lambda x, y: (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def det(M, rows, cols):
    total = ZERO
    for p in itertools.permutations(range(len(cols))):
        term = (-1 if sum(p[j] > p[i] for i in range(len(p)) for j in range(i)) % 2 else 1, 0)
        for i, j in enumerate(p):
            term = mul(term, M[rows[i]][cols[j]])
        total = add(total, term)
    return total


def matmul(X, Y):
    n = len(X)
    out = [[ZERO] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                out[i][j] = add(out[i][j], mul(X[i][k], Y[k][j]))
    return out


def gate(n, l, r, a, b, c, d, diag=None):
    M = [[(diag[i] if diag else ONE) if i == j else ZERO for j in range(n)] for i in range(n)]
    M[l][l], M[l][r], M[r][l], M[r][r] = a, b, c, d
    return M


entry = lambda: (rng.randint(-2, 2), rng.randint(-2, 2))
nonzero = lambda: next(z for z in iter(entry, None) if z != ZERO)
subsets = lambda n, k: [tuple(s) for s in itertools.combinations(range(n), k)]
swap = lambda S, l, r: tuple(sorted({r if x == l else l if x == r else x for x in S}))

# A and B
pairs = nz = nz_swap = nz_coll = sharp = 0
for n in range(2, 6):
    for l, r in itertools.permutations(range(n), 2):
        for trial in range(6):
            diag = [entry() for _ in range(n)] if trial % 2 else None
            M = gate(n, l, r, entry(), entry(), entry(), entry(), diag)
            for k in range(n + 1):
                for S in subsets(n, k):
                    for T in subsets(n, k):
                        pairs += 1
                        if det(M, T, S) != ZERO:
                            nz += 1
                            assert T == S or T == swap(S, l, r), (n, l, r, S, T)
                            nz_swap += T != S
                            nz_coll += (l in S) == (r in S)
            while True:
                a, b, c, d = nonzero(), nonzero(), nonzero(), nonzero()
                if add(mul(a, d), mul((-1, 0), mul(b, c))) != ZERO:
                    break
            G = gate(n, l, r, a, b, c, d)
            for k in range(n + 1):
                for S in subsets(n, k):
                    for T in {S, swap(S, l, r)}:
                        assert det(G, T, S) != ZERO
                        sharp += 1
print(f'A: {pairs} square minors, {nz} nonzero ({nz_swap} swaps, {nz_coll} in collision sectors): all stay or swap')
print(f'B: {sharp} legal stay/swap minors of generic gates: all nonzero')


def history_counts(n, sched, k):
    """For each (S0, Sm): number of nonzero histories and the sum of their amplitudes."""
    result = {}
    for S0 in subsets(n, k):
        layer = {S0: (1, ONE)}
        for (l, r), G in sched:
            nxt = {}
            for h, (cnt, amp) in layer.items():
                for U in subsets(n, k):
                    f = det(G, U, h)
                    if f != ZERO:
                        c0, a0 = nxt.get(U, (0, ZERO))
                        nxt[U] = (c0 + cnt, add(a0, mul(amp, f)))
            layer = nxt
        for Sm, val in layer.items():
            result[(S0, Sm)] = val
    return result


def run_schedule(n, edges):
    sched = [((l, r), gate(n, l, r, entry(), entry(), entry(), entry())) for l, r in edges]
    P = gate(n, 0, 1, ONE, ZERO, ZERO, ONE)
    for _, G in sched:
        P = matmul(G, P)
    worst = 0
    for k in range(n + 1):
        for (S0, Sm), (cnt, amp) in history_counts(n, sched, k).items():
            assert amp == det(P, Sm, S0), 'Cauchy-Binet mismatch'
            worst = max(worst, cnt)
    return worst


# C
schedules = 0
for n in range(2, 6):
    for trial in range(12):
        order = list(range(n)); rng.shuffle(order)
        edges = [(order[i], rng.choice(order[:i])) for i in range(1, n) if rng.random() < 0.85]
        rng.shuffle(edges)
        edges = [(r, l) if rng.random() < 0.5 else (l, r) for l, r in edges]
        assert run_schedule(n, edges) <= 1, edges
        schedules += 1
print(f'C: {schedules} random forest schedules (n=2..5): at most one nonzero history per endpoint pair; Cauchy-Binet sums match')

# D
tri = max(run_schedule(3, [(0, 1), (1, 2), (2, 0)]) for _ in range(5))
rep = max(run_schedule(2, [(0, 1), (0, 1)]) for _ in range(5))
assert tri >= 2 and rep >= 2
print(f'D: controls without the bridge hypothesis reach {tri} (triangle) and {rep} (repeated edge) nonzero histories')

if '--lean' not in sys.argv:
    sys.exit(0)

PROBE = r'''import ForestGateSupport
open MathCommons.ForestTrajectory MathCommons.ForestGateSupport
def verts : List (Fin 3) := [0, 1, 2]
def e01 : Edge verts where
  nodup := by decide
  complete := by decide
  left := 0
  right := 1
  distinct := by decide
  rest := [2]
  enumerate := List.Perm.refl _
  avoids := by decide
def e12 : Edge verts where
  nodup := by decide
  complete := by decide
  left := 1
  right := 2
  distinct := by decide
  rest := [0]
  enumerate := by decide
  avoids := by decide
def pathEdges : Nat → Edge verts
  | 1 => e12
  | _ => e01
theorem side0 (x y : Fin 3) (h : ReachWithout pathEdges 2 0 x y) : (x = 0 ↔ y = 0) := by
  induction h with
  | refl => exact Iff.rfl
  | edge k within different =>
    have : k = 1 := by omega
    subst this; decide
  | symm _ ih => exact ih.symm
  | trans _ _ ih1 ih2 => exact ih1.trans ih2
theorem side1 (x y : Fin 3) (h : ReachWithout pathEdges 2 1 x y) : (x = 2 ↔ y = 2) := by
  induction h with
  | refl => exact Iff.rfl
  | edge k within different =>
    have : k = 0 := by omega
    subst this; decide
  | symm _ ih => exact ih.symm
  | trans _ _ ih1 ih2 => exact ih1.trans ih2
theorem pathBridges : EveryEdgeIsBridge pathEdges 2 := by
  intro j hj reach
  match j, hj with
  | 0, _ => exact absurd ((side0 _ _ reach).mp (by decide) : (pathEdges 0).right = 0) (by decide)
  | 1, _ => exact absurd ((side1 _ _ reach).mpr (by decide) : (pathEdges 1).left = 2) (by decide)
def swapGate := twoCoordinateGate gaussianArithmetic e01 (0, 0) (1, 0) (1, 0) (0, 0)
def sIn : Subset (Fin 3) := fun x => x == 0
def tOut : Subset (Fin 3) := fun x => x == 1
theorem swapMinorNonzero :
    minorDet gaussianArithmetic swapGate (verts.filter tOut) (verts.filter sIn) ≠ (0, 0) := by
  have h : minorDet gaussianArithmetic swapGate (verts.filter tOut) (verts.filter sIn) = (1, 0) := by
    simp [minorDet, columnChoices, scalarSum, signed, swapGate, twoCoordinateGate,
      gaussianArithmetic, verts, sIn, tOut, e01]
    exact ⟨by rw [Rat.sub_eq_add_neg]; show (1 : Rat) + 0 + 0 = 1; rw [Rat.add_zero, Rat.add_zero],
      by rw [Rat.add_zero, Rat.add_zero]⟩
  rw [h]; decide
theorem swapStep : Step e01 sIn tOut := gaussianGateMinorStep e01 _ _ _ _ sIn tOut swapMinorNonzero
theorem swapIsExchange : tOut = swapSubset e01 sIn ∧ tOut ≠ sIn :=
  ⟨by funext x; revert x; decide, fun h => by have := congrFun h 1; revert this; decide⟩
#print axioms pathBridges
#print axioms swapStep
#print axioms swapIsExchange
def matA : Nat → Nat → Rat := fun i j => ((MATRIX : List (List Int)).getD i []).getD j 0
EVALS
'''
A = [[3, -1, 4, 1, -5], [9, 2, -6, 5, 3], [-5, 8, 9, -7, 9], [3, 2, -3, 8, 4], [-6, 2, 6, 4, -3]]
MINORS = [([0, 1, 2, 3, 4], [0, 1, 2, 3, 4]), ([0, 2, 3], [1, 2, 4]), ([4, 1], [3, 0]),
          ([1, 2, 3, 4], [0, 1, 3, 4]), ([0, 1], [0, 1, 2]), ([0, 1, 2], [2, 3]), ([], [])]
Ag = [[(x, 0) for x in row] for row in A]
expected = [det(Ag, r, c)[0] if len(r) == len(c) else 0 for r, c in MINORS]
probe = PROBE.replace('MATRIX', str(A)).replace('EVALS', '\n'.join(
    f'#eval minorDet rationalArithmetic matA {r} {c}' for r, c in MINORS))
formal = Path(__file__).resolve().parents[3] / 'formal'
with tempfile.TemporaryDirectory() as build:
    env = {**os.environ, 'LEAN_PATH': build}
    for mod in ('ForestCut', 'ForestTrajectory', 'ForestGateSupport'):
        subprocess.run(['lean', '-DwarningAsError=true', '-o', f'{build}/{mod}.olean', f'{mod}.lean'],
                       cwd=formal, env=env, check=True, capture_output=True, timeout=300)
    Path(build, 'Probe.lean').write_text(probe)
    run = subprocess.run(['lean', '-DwarningAsError=true', 'Probe.lean'], cwd=build, env=env,
                         capture_output=True, text=True, timeout=300)
assert run.returncode == 0, run.stdout + run.stderr
axioms = re.findall(r"depends on axioms: \[([^\]]*)\]", run.stdout)
assert len(axioms) == 3 and all(set(a.split(', ')) <= {'propext', 'Classical.choice', 'Quot.sound'} for a in axioms)
got = [int(x) for x in run.stdout.split('\n') if re.fullmatch(r'-?\d+', x.strip())]
assert got == expected, (got, expected)
print(f'E: Lean probe compiled; bridge and swap examples use only standard axioms; minorDet = Leibniz on {len(MINORS)} minors {got}')
