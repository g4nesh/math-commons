"""Kernel-replay Lemma L, enforce its exact type, and audit every public axiom set."""
from pathlib import Path
from tempfile import TemporaryDirectory
import hashlib
import json
import os
import re
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOURCE = ROOT / 'formal/OrthogonalCube.lean'
ALLOWED = {'propext', 'Classical.choice', 'Quot.sound'}
DECLARATIONS = [
    'squareNonnegative', 'quadraticNonnegative', 'lagrangeIdentity',
    'quadraticCauchy', 'quadraticReverseCauchy', 'firstPartNonnegative',
    'firstPartCauchy', 'secondPartReverseCauchy', 'complementPolarization',
    'totalAtOne', 'strictOrderOfSquares', 'weightedCauchy',
    'strictWeightedCross', 'contradictionFromParts', 'lemmaLOrdered',
    'cubeFormSwap', 'lemmaL', 'lemmaLMax',
]
PREFIX = 'MathCommons.OrthogonalCube.'
EXPECTED = {PREFIX + name for name in DECLARATIONS}
TYPE_CHECK = '''import OrthogonalCube
open Std Lean.Grind MathCommons.OrthogonalCube
example {R : Type} [Field R] [LE R] [LT R]
    [LawfulOrderLT R] [IsLinearOrder R] [OrderedRing R]
    (p q r u a b c d : R)
    (hp : 0 ≤ p) (hq : 0 ≤ q) (hr : 0 ≤ r) (hu : 0 ≤ u)
    (balance : p+q=r+u)
    (ha : 0 ≤ a ∧ a ≤ 1) (hb : 0 ≤ b ∧ b ≤ 1)
    (hc : 0 ≤ c ∧ c ≤ 1) (hd : 0 ≤ d ∧ d ≤ 1) :
    0 ≤ p*r*(a*a+c*c)+q*u*(b*b+d*d)-2*p*q*a*b-2*r*u*c*d ∨
    0 ≤ p*r*((1-a)*(1-a)+(1-c)*(1-c))+q*u*((1-b)*(1-b)+(1-d)*(1-d))
      -2*p*q*(1-a)*(1-b)-2*r*u*(1-c)*(1-d) :=
  lemmaL p q r u a b c d hp hq hr hu balance ha hb hc hd
example (p q r u a b c d : Rat)
    (hp : 0 ≤ p) (hq : 0 ≤ q) (hr : 0 ≤ r) (hu : 0 ≤ u)
    (balance : p+q=r+u)
    (ha : 0 ≤ a ∧ a ≤ 1) (hb : 0 ≤ b ∧ b ≤ 1)
    (hc : 0 ≤ c ∧ c ≤ 1) (hd : 0 ≤ d ∧ d ≤ 1) :
    0 ≤ orderedMax (cubeForm p q r u a b c d)
      (cubeForm p q r u (1-a) (1-b) (1-c) (1-d)) :=
  lemmaLMax p q r u a b c d hp hq hr hu balance ha hb hc hd
'''


def check_reports(output):
    matches = re.findall(r"'([\w.]+)' (?:does not depend on any axioms|depends on axioms:\s*\[([^\]]*)\])", output)
    seen = set()
    for name, deps in matches:
        used = {part.strip() for part in deps.split(',') if part.strip()}
        if used - ALLOWED:
            raise SystemExit(f'{name}: unapproved axioms {sorted(used - ALLOWED)}')
        if name in seen:
            raise SystemExit(f'Duplicate axiom report: {name}')
        seen.add(name)
    if seen != EXPECTED:
        raise SystemExit(f'Axiom report mismatch: missing={sorted(EXPECTED-seen)}, extra={sorted(seen-EXPECTED)}')


def main():
    lean = sys.argv[1] if len(sys.argv) == 2 else 'lean'
    version = (ROOT / 'formal/lean-toolchain').read_text().strip().split(':v')[1]
    reported = subprocess.check_output([lean, '--version'], text=True)
    if not re.search(r'\bversion ' + re.escape(version) + r'[, )]', reported):
        raise SystemExit(f'Expected Lean {version}, found {reported.strip()}')
    audit = json.loads((HERE / 'audit.json').read_text())
    digest = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    if audit['sourceSha256'] != digest or audit['declarations'] != DECLARATIONS:
        raise SystemExit('Source hash or declaration inventory differs from audit.json')
    if audit['leanVersion'] != version or set(audit['allowedAxioms']) != ALLOWED:
        raise SystemExit('Audit version or axiom policy differs from checker')
    with TemporaryDirectory(prefix='commons-orthogonal-cube-') as temp:
        env = {**os.environ, 'LEAN_PATH': temp}
        run = subprocess.run([lean, '-DwarningAsError=true', '-o', str(Path(temp, 'OrthogonalCube.olean')), str(SOURCE)],
                             cwd=ROOT, env=env, text=True, capture_output=True, timeout=120)
        print(run.stdout, end='')
        if run.returncode:
            raise SystemExit(run.stderr or f'Lean exited {run.returncode}')
        check_reports(run.stdout)
        type_file = Path(temp, 'StatementCheck.lean')
        type_file.write_text(TYPE_CHECK)
        typed = subprocess.run([lean, '-DwarningAsError=true', str(type_file)],
                               cwd=temp, env=env, text=True, capture_output=True, timeout=60)
        if typed.returncode:
            raise SystemExit(typed.stdout + typed.stderr)
    print(f'PASS: {len(EXPECTED)} kernel-checked declarations; exact generic statement and Rat instance.')
    print(f'SHA256 {digest}')


if __name__ == '__main__':
    main()
