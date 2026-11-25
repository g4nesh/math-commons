"""Independent review checks for task mdo-cut-subset-bridge (formal/ForestTrajectory.lean).

Default run (Python 3 with numpy, a few seconds):
  [A] gate support: for random two-coordinate unitary gates in dimensions 2..5, every
      nonzero minor det G[R', R] has R' = R or R' = R with the two endpoints exchanged,
      and the squared magnitudes follow the baselines transition table (seed 20261125);
  [B] exhaustive model check: for every edge sequence on d <= 4 vertices with m <= 5
      edges, the Lean bridge hypothesis (deleting edge j disconnects its endpoints)
      equals acyclicity, and under it every (start, end) pair has at most one legal
      stay-or-swap subset history.
With --lean (needs `lean` 4.24.0 on PATH): checks the SHA-256 of the two pinned Lean
files, compiles them with warnings as errors in a temporary directory, checks the 17
axiom reports, and compiles the reviewer's own Lean audit below against them.
"""
import hashlib
import itertools
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

PINNED = {
    'ForestCut.lean': '0bfc29b6c4470e3c11cfdbdb561616f55bc3d73af78d53e8edc68feb44b1b5ae',
    'ForestTrajectory.lean': 'b5dcffb2a165702ef46a2a67e5a7d9049a5f5ba3af3cfd9cc2491af5e0fc5910',
}
ALLOWED = {'propext', 'Classical.choice', 'Quot.sound'}


def gate_support(seed=20261125, trials=4):
    rng = np.random.default_rng(seed)
    checked = 0
    for d in range(2, 6):
        for p, q in itertools.combinations(range(d), 2):
            for _ in range(trials):
                z = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
                g, _ = np.linalg.qr(z)
                G = np.eye(d, dtype=complex)
                G[np.ix_([p, q], [p, q])] = g
                qr = abs(g[0, 1]) ** 2
                for k in range(d + 1):
                    for R in itertools.combinations(range(d), k):
                        swap = tuple(sorted({p: q, q: p}.get(x, x) for x in R))
                        for Rp in itertools.combinations(range(d), k):
                            m = np.linalg.det(G[np.ix_(Rp, R)]) if k else 1.0
                            hit = len({p, q} & set(R))
                            if Rp == R:
                                want = 1 - qr if hit == 1 else 1.0
                            elif Rp == swap:
                                want = qr
                            else:
                                want = 0.0
                            assert abs(abs(m) ** 2 - want) < 1e-12, (d, p, q, R, Rp)
                            checked += 1
    print(f'[A] gate support and transition table: {checked} minors checked, all as stated')


def model_check(max_d=4, max_m=5):
    def find(par, x):
        while par[x] != x:
            x = par[x]
        return x

    def connected(d, es, u, v):
        par = list(range(d))
        for a, b in es:
            par[find(par, a)] = find(par, b)
        return find(par, u) == find(par, v)

    total = bridge_count = nonunique_cyclic = 0
    for d in range(2, max_d + 1):
        pairs = list(itertools.combinations(range(d), 2))
        for m in range(max_m + 1):
            for es in itertools.product(pairs, repeat=m):
                bridge = all(not connected(d, es[:j] + es[j + 1:], *es[j]) for j in range(m))
                forest = all(not connected(d, es[:j], *es[j]) for j in range(m))
                assert bridge == forest, (d, es)
                count = {(s, s): 1 for s in range(1 << d)}
                for a, b in es:
                    nxt = {}
                    for (s, cur), c in count.items():
                        sw = cur
                        if ((cur >> a) & 1) != ((cur >> b) & 1):
                            sw = cur ^ (1 << a) ^ (1 << b)
                        for t in {cur, sw}:
                            nxt[(s, t)] = nxt.get((s, t), 0) + c
                    count = nxt
                unique = max(count.values()) == 1
                assert unique or not bridge, (d, es)
                total += 1
                bridge_count += bridge
                nonunique_cyclic += (not bridge) and (not unique)
    print(f'[B] {total} edge sequences (d <= {max_d}, m <= {max_m}): bridge == forest in all; '
          f'{bridge_count} forests, all with unique histories; '
          f'{nonunique_cyclic} of {total - bridge_count} non-forests have two histories with equal endpoints')


def lean_check(root):
    formal = root / 'formal'
    for name, want in PINNED.items():
        got = hashlib.sha256((formal / name).read_bytes()).hexdigest()
        assert got == want, f'{name}: sha256 {got} differs from the reviewed file'
    with tempfile.TemporaryDirectory(prefix='review-cut-subset-') as tmp:
        env = {**os.environ, 'LEAN_PATH': tmp}
        reports = {}
        for name in ('ForestCut.lean', 'ForestTrajectory.lean'):
            run = subprocess.run(['lean', '-DwarningAsError=true', '-o',
                                  str(Path(tmp, name[:-5] + '.olean')), name],
                                 cwd=formal, env=env, text=True, capture_output=True, timeout=300)
            assert run.returncode == 0, run.stdout + run.stderr
            reports[name] = re.findall(r"'([\w.]+)' depends on axioms: \[([^\]]*)\]", run.stdout)
        traj = reports['ForestTrajectory.lean']
        assert len(traj) == 17, len(traj)
        for name, deps in sum(reports.values(), []):
            assert {x.strip() for x in deps.split(',')} <= ALLOWED, (name, deps)
        audit = Path(tmp, 'ReviewAudit.lean')
        audit.write_text(AUDIT, encoding='utf-8')
        run = subprocess.run(['lean', '-DwarningAsError=true', str(audit)], cwd=tmp, env=env,
                             text=True, capture_output=True, timeout=300)
        assert run.returncode == 0, run.stdout + run.stderr
        own = re.findall(r"'([\w.]+)' depends on axioms: \[([^\]]*)\]", run.stdout)
        assert len(own) == 4 and all({x.strip() for x in d.split(',')} <= ALLOWED for _, d in own), own
    print('[L] sha256 pinned; ForestCut and ForestTrajectory compile with warnings as errors; '
          f'{len(traj)} ForestTrajectory axiom reports within {sorted(ALLOWED)}; '
          'reviewer audit (restatement, path instance, triangle control) compiles')


AUDIT = r'''import ForestTrajectory
open MathCommons.ForestTrajectory

namespace ReviewAudit

/-! (1) Statement pin: the theorem restated with a textbook walk notion and the
stay-or-swap rule written out, so it depends on `Step`, `swapSubset`, `Trajectory`
and `ReachWithout` only through definitional unfolding and the lemma `toWalk`. -/

/-- A walk from `x` to `z` along indexed edges `k < n`, `k ≠ j`, either direction. -/
inductive Walk {V : Type} {vs : List V} (edges : Nat → Edge vs) (n j : Nat) : V → V → Prop
  | nil (x : V) : Walk edges n j x x
  | cons {x y z : V} (k : Nat) (hk : k < n) (hne : k ≠ j)
      (adj : (x = (edges k).left ∧ y = (edges k).right) ∨
             (x = (edges k).right ∧ y = (edges k).left))
      (rest : Walk edges n j y z) : Walk edges n j x z

theorem Walk.append {V : Type} {vs : List V} {edges : Nat → Edge vs} {n j : Nat}
    {x y z : V} (p : Walk edges n j x y) (q : Walk edges n j y z) : Walk edges n j x z := by
  induction p with
  | nil => exact q
  | cons k hk hne adj _ ih => exact .cons k hk hne adj (ih q)

theorem Walk.reverse {V : Type} {vs : List V} {edges : Nat → Edge vs} {n j : Nat}
    {x y : V} (p : Walk edges n j x y) : Walk edges n j y x := by
  induction p with
  | nil x => exact .nil x
  | cons k hk hne adj _ ih =>
    refine ih.append (.cons k hk hne ?_ (.nil _))
    rcases adj with ⟨h1, h2⟩ | ⟨h1, h2⟩
    · exact Or.inr ⟨h2, h1⟩
    · exact Or.inl ⟨h2, h1⟩

theorem toWalk {V : Type} {vs : List V} {edges : Nat → Edge vs} {n j : Nat} {x y : V}
    (r : ReachWithout edges n j x y) : Walk edges n j x y := by
  induction r with
  | refl x => exact .nil x
  | edge k hk hne => exact .cons k hk hne (Or.inl ⟨rfl, rfl⟩) (.nil _)
  | symm _ ih => exact ih.reverse
  | trans _ _ ih1 ih2 => exact ih1.append ih2

theorem restated {V : Type} [DecidableEq V] {vs : List V} (edges : Nat → Edge vs) (n : Nat)
    (bridge : ∀ j, j < n → ¬ Walk edges n j (edges j).left (edges j).right)
    (a b : Nat → V → Bool)
    (ha : ∀ k, k < n → a (k + 1) = a k ∨ a (k + 1) = fun x =>
      if x = (edges k).left then a k (edges k).right
      else if x = (edges k).right then a k (edges k).left else a k x)
    (hb : ∀ k, k < n → b (k + 1) = b k ∨ b (k + 1) = fun x =>
      if x = (edges k).left then b k (edges k).right
      else if x = (edges k).right then b k (edges k).left else b k x)
    (h0 : a 0 = b 0) (hn : a n = b n) : ∀ k, k ≤ n → a k = b k :=
  forestTrajectoryUnique edges n (fun j hj r => bridge j hj (toWalk r)) a b ha hb h0 hn

/-! (2) Non-vacuity on the path 0-1-2, and (3) a triangle control. -/

def vs : List (Fin 3) := [0, 1, 2]

def mk (l r o : Fin 3) (hd : l ≠ r) (hp : vs.Perm [l, r, o])
    (ha : ∀ x ∈ [o], x ≠ l ∧ x ≠ r) : Edge vs :=
  { nodup := by decide, complete := by decide, left := l, right := r, distinct := hd,
    rest := [o], enumerate := hp, avoids := ha }

def e01 : Edge vs := mk 0 1 2 (by decide) (by decide) (by decide)
def e12 : Edge vs := mk 1 2 0 (by decide) (by decide) (by decide)
def e02 : Edge vs := mk 0 2 1 (by decide) (by decide) (by decide)

def path : Nat → Edge vs | 0 => e01 | _ => e12
def tri : Nat → Edge vs | 0 => e01 | 1 => e12 | _ => e02

theorem pathBridges : EveryEdgeIsBridge path 2 := by
  intro j hj
  match j, hj with
  | 0, _ =>
    intro r
    have inv : ∀ x y, ReachWithout path 2 0 x y → (x = 0 ↔ y = 0) := by
      intro x y h
      induction h with
      | refl => exact Iff.rfl
      | edge k hk hne =>
        obtain rfl : k = 1 := by omega
        decide
      | symm _ ih => exact ih.symm
      | trans _ _ ih1 ih2 => exact ih1.trans ih2
    exact absurd (inv _ _ r) (by decide)
  | 1, _ =>
    intro r
    have inv : ∀ x y, ReachWithout path 2 1 x y → (x = 2 ↔ y = 2) := by
      intro x y h
      induction h with
      | refl => exact Iff.rfl
      | edge k hk hne =>
        obtain rfl : k = 0 := by omega
        decide
      | symm _ ih => exact ih.symm
      | trans _ _ ih1 ih2 => exact ih1.trans ih2
    exact absurd (inv _ _ r) (by decide)

def single (v : Fin 3) : Subset (Fin 3) := fun x => decide (x = v)

/-- Any legal history on the path from {0} to {2} must pass through {1}. -/
theorem pathForcesMiddle (b : Nat → Subset (Fin 3)) (hb : Trajectory path 2 b)
    (h0 : b 0 = single 0) (h2 : b 2 = single 2) : b 1 = single 1 := by
  let a : Nat → Subset (Fin 3) := fun k => match k with | 0 => single 0 | 1 => single 1 | _ => single 2
  have ha : Trajectory path 2 a := by
    intro k hk
    match k, hk with
    | 0, _ => exact Or.inr (by funext x; revert x; decide)
    | 1, _ => exact Or.inr (by funext x; revert x; decide)
  exact (forestTrajectoryUnique path 2 pathBridges a b ha hb h0.symm h2.symm 1 (by decide)).symm

/-- The triangle violates the bridge hypothesis, and two legal histories with equal
endpoints differ at time 1, so the hypothesis cannot be dropped. -/
theorem triangleControl : ¬ EveryEdgeIsBridge tri 3 ∧
    ∃ a b : Nat → Subset (Fin 3), Trajectory tri 3 a ∧ Trajectory tri 3 b ∧
      a 0 = b 0 ∧ a 3 = b 3 ∧ a 1 ≠ b 1 := by
  refine ⟨fun h => h 0 (by decide) ?_, ?_⟩
  · exact .trans (.edge 2 (by decide) (by decide)) (.symm (.edge 1 (by decide) (by decide)))
  · refine ⟨fun k => match k with | 0 => single 0 | 1 => single 1 | _ => single 2,
            fun k => match k with | 0 => single 0 | 1 => single 0 | 2 => single 0 | _ => single 2,
            ?_, ?_, rfl, rfl, ?_⟩
    · intro k hk
      match k, hk with
      | 0, _ => exact Or.inr (by funext x; revert x; decide)
      | 1, _ => exact Or.inr (by funext x; revert x; decide)
      | 2, _ => exact Or.inl rfl
    · intro k hk
      match k, hk with
      | 0, _ => exact Or.inl rfl
      | 1, _ => exact Or.inl rfl
      | 2, _ => exact Or.inr (by funext x; revert x; decide)
    · intro h; have := congrFun h 1; revert this; decide

#check @forestTrajectoryUnique
#print axioms restated
#print axioms pathBridges
#print axioms pathForcesMiddle
#print axioms triangleControl

end ReviewAudit
'''

if __name__ == '__main__':
    args = sys.argv[1:]
    root = Path(__file__).resolve().parents[3]
    if '--root' in args:
        root = Path(args[args.index('--root') + 1]).resolve()
    gate_support()
    model_check()
    if '--lean' in args:
        lean_check(root)
