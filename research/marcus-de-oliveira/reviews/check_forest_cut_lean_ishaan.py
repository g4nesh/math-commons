#!/usr/bin/env python3
"""Independent checks of formal/ForestCut.lean (checkpoint mdo-forest-cut-lean-published-20261125).

Review aid for mdo-forest-cut-lean, written by ishaan-air. Standard library only.

1. The reviewed file has SHA-256 0bfc29b6...b1b5ae, the hash cited by the checkpoint.
2. Lean 4.24.0 accepts it with warnings as errors, and the eight published theorems
   report only the axioms propext and Quot.sound.
3. Lean accepts the file followed by a reviewer block that (a) pins every definition and
   theorem statement to an independent restatement, (b) instantiates all hypotheses on a
   concrete history with a real jump, (c) shows the one-opportunity hypothesis is needed,
   and (d) notes that zero transitions give no opportunity.
4. Bridge test, outside Lean: on fixed and seeded random forests (n <= 7 vertices, any
   gate order), every stay/exchange subset trajectory gives, for every edge and both sides
   of its cut, a count history satisfying SingleOpportunity with jump in {-1,0,1}, and no
   two trajectories share both endpoints. A triangle and a repeated edge fail uniqueness.
   This is finite evidence about the unformalized bridge, not a proof of it.

Usage, from the repository root:  python3 -I research/marcus-de-oliveira/reviews/check_forest_cut_lean_ishaan.py [lean]
"""
import hashlib, random, re, subprocess, sys, tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SOURCE = ROOT / 'formal' / 'ForestCut.lean'
SHA = '0bfc29b6c4470e3c11cfdbdb561616f55bc3d73af78d53e8edc68feb44b1b5ae'
NAMES = ['telescoping', 'constantOnInterval', 'singleCutChange', 'totalChangeAtOpportunity',
         'sameEndpointsSameCrossing', 'cutHistoryUnique', 'crossingDeltaInjective', 'crossingDecisionUnique']
ALLOWED = {'propext', 'Quot.sound'}
LEAN = sys.argv[1] if len(sys.argv) > 1 else 'lean'

BLOCK = r'''
namespace ReviewForestCut
open MathCommons.ForestCut
example (f : Nat → Int) : totalChange f 0 = 0 := rfl
example (f : Nat → Int) (n : Nat) : totalChange f (n + 1) = totalChange f n + (f (n + 1) - f n) := rfl
example (f : Nat → Int) (n j : Nat) :
    SingleOpportunity f n j ↔ (j < n ∧ ∀ k, k < n → k ≠ j → f (k + 1) = f k) := Iff.rfl
example : crossingDelta .stay = 0 ∧ crossingDelta .enter = 1 ∧ crossingDelta .leave = -1 := ⟨rfl, rfl, rfl⟩
theorem pinTelescoping : ∀ (f : Nat → Int) (n : Nat), totalChange f n = f n - f 0 := telescoping
theorem pinConstant : ∀ {α : Type} (f : Nat → α) (a b : Nat), a ≤ b →
    (∀ k, a ≤ k → k < b → f (k + 1) = f k) → f b = f a := @constantOnInterval
theorem pinSingle : ∀ (f : Nat → Int) (n j : Nat), SingleOpportunity f n j →
    f (j + 1) - f j = f n - f 0 := singleCutChange
theorem pinTotal : ∀ (f : Nat → Int) (n j : Nat), SingleOpportunity f n j →
    totalChange f n = f (j + 1) - f j := totalChangeAtOpportunity
theorem pinSame : ∀ (f g : Nat → Int) (n j : Nat), SingleOpportunity f n j → SingleOpportunity g n j →
    f 0 = g 0 → f n = g n → f (j + 1) - f j = g (j + 1) - g j := sameEndpointsSameCrossing
theorem pinHistory : ∀ (f g : Nat → Int) (n j : Nat), SingleOpportunity f n j → SingleOpportunity g n j →
    f 0 = g 0 → f n = g n → ∀ k, k ≤ n → f k = g k := cutHistoryUnique
theorem pinInjective : ∀ x y : Crossing, crossingDelta x = crossingDelta y → x = y := crossingDeltaInjective
theorem pinDecision : ∀ (f g : Nat → Int) (n j : Nat), SingleOpportunity f n j → SingleOpportunity g n j →
    f 0 = g 0 → f n = g n → ∀ x y : Crossing,
    f (j + 1) - f j = crossingDelta x → g (j + 1) - g j = crossingDelta y → x = y := crossingDecisionUnique
def hEnter (k : Nat) : Int := if k ≤ 2 then 2 else 3
def gEnter (k : Nat) : Int := if k < 3 then 2 else 3
theorem hOpp : SingleOpportunity hEnter 5 2 := by unfold SingleOpportunity; decide
theorem gOpp : SingleOpportunity gEnter 5 2 := by unfold SingleOpportunity; decide
theorem lastOpp : SingleOpportunity hEnter 3 2 := by unfold SingleOpportunity; decide
theorem hJump : hEnter 3 - hEnter 2 = 1 := singleCutChange hEnter 5 2 hOpp
theorem hgAgree : ∀ k, k ≤ 5 → hEnter k = gEnter k :=
  cutHistoryUnique hEnter gEnter 5 2 hOpp gOpp (by decide) (by decide)
theorem decisionIsEnter (y : Crossing) (hy : gEnter 3 - gEnter 2 = crossingDelta y) : y = .enter :=
  (crossingDecisionUnique hEnter gEnter 5 2 hOpp gOpp (by decide) (by decide) .enter y (by decide) hy).symm
def bump (k : Nat) : Int := if k = 1 then 1 else 0
def flat (_ : Nat) : Int := 0
theorem hypothesisNeeded : bump 0 = flat 0 ∧ bump 2 = flat 2 ∧ bump 1 - bump 0 ≠ flat 1 - flat 0 ∧
    ¬ SingleOpportunity bump 2 0 ∧ ¬ SingleOpportunity bump 2 1 := by unfold SingleOpportunity; decide
theorem noOpportunityAtZero (f : Nat → Int) (j : Nat) : ¬ SingleOpportunity f 0 j := fun h => Nat.not_lt_zero j h.1
#print axioms hgAgree
#print axioms decisionIsEnter
#print axioms hypothesisNeeded
end ReviewForestCut
'''


def run_lean(text):
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / 'Check.lean'
        path.write_text(text, encoding='utf-8')
        run = subprocess.run([LEAN, '-DwarningAsError=true', path.name], cwd=tmp, text=True,
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=300)
    if run.returncode:
        sys.exit('Lean failed:\n' + run.stdout)
    return {name: {a.strip() for a in deps.split(',') if a.strip()} for name, deps in
            re.findall(r"'([\w.]+)' (?:does not depend on any axioms|depends on axioms: \[([^\]]*)\])", run.stdout)}


def trajectories(n, edges, start):
    paths = [[start]]
    for i, j in edges:
        nxt = []
        for p in paths:
            x = p[-1]
            nxt.append(p + [x])
            if ((x >> i) ^ (x >> j)) & 1:  # exactly one endpoint occupied: exchange allowed
                nxt.append(p + [x ^ (1 << i) ^ (1 << j)])
        paths = nxt
    return paths


def side(n, edges, r):
    """Vertices reachable from the first endpoint of edge r without using edge r."""
    seen, stack = {edges[r][0]}, [edges[r][0]]
    while stack:
        v = stack.pop()
        for s, (a, b) in enumerate(edges):
            for u, w in ((a, b), (b, a)):
                if s != r and u == v and w not in seen:
                    seen.add(w)
                    stack.append(w)
    return seen


def bridge_check(n, edges):
    """Return (number of endpoint pairs with >1 trajectory, number of cut histories checked)."""
    clashes, histories = 0, 0
    sides = [side(n, edges, r) for r in range(len(edges))]
    for t in range(1 << n):
        ends = {}
        for p in trajectories(n, edges, t):
            ends[p[-1]] = ends.get(p[-1], 0) + 1
            for r, c in enumerate(sides):
                if edges[r][1] in c:
                    continue  # edge r lies on a cycle or is repeated: no cut
                for cut in (c, set(range(n)) - c):
                    f = [sum((x >> v) & 1 for v in cut) for x in p]
                    assert all(f[k + 1] == f[k] for k in range(len(edges)) if k != r), (edges, r, p)
                    assert f[r + 1] - f[r] in (-1, 0, 1)
                    histories += 1
        clashes += sum(c > 1 for c in ends.values())
    return clashes, histories


def random_forest(rng, n):
    edges = [(rng.randrange(v), v) for v in range(1, n) if rng.random() < 0.8]
    edges = [(b, a) if rng.random() < 0.5 else (a, b) for a, b in edges]
    rng.shuffle(edges)
    return edges


def main():
    data = SOURCE.read_bytes()
    if hashlib.sha256(data).hexdigest() != SHA:
        sys.exit(f'{SOURCE} is not the reviewed file')
    version = subprocess.check_output([LEAN, '--version'], text=True)
    if 'version 4.24.0' not in version:
        sys.exit('Expected Lean 4.24.0, found ' + version.strip())
    text = data.decode('utf-8')
    for label, src in (('file', text), ('file+review', text + BLOCK)):
        reports = run_lean(src)
        missing = [n for n in NAMES if 'MathCommons.ForestCut.' + n not in reports]
        extra = {n: u - ALLOWED for n, u in reports.items() if u - ALLOWED}
        if missing or extra:
            sys.exit(f'{label}: missing reports {missing}, unapproved axioms {extra}')
        print(f'Lean {label}: accepted, {len(reports)} axiom reports, all within {sorted(ALLOWED)}')
    rng = random.Random(20261125)
    forests = [(1, []), (6, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)]), (6, [(4, 5), (2, 3), (0, 1), (1, 2), (3, 4)]),
               (6, [(0, 5), (0, 1), (0, 4), (0, 2), (0, 3)])]
    forests += [(n, random_forest(rng, n)) for n in rng.choices(range(2, 8), k=300)]
    total = 0
    for n, edges in forests:
        clashes, histories = bridge_check(n, edges)
        assert clashes == 0, (n, edges)
        total += histories
    print(f'Bridge: {len(forests)} forests, {total} cut histories satisfy SingleOpportunity, no endpoint clash')
    for name, n, edges in (('triangle', 3, [(0, 1), (1, 2), (0, 2)]), ('repeated edge', 2, [(0, 1), (0, 1)])):
        clashes, _ = bridge_check(n, edges)
        assert clashes > 0
        print(f'Control {name}: {clashes} endpoint pairs with two trajectories (hypothesis fails, as expected)')
    print('All checks passed.')


if __name__ == '__main__':
    main()
