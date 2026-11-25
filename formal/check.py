"""Recheck all listed public theorems and their transitive axiom reports."""
from pathlib import Path
import re
import subprocess
import sys

root = Path(__file__).resolve().parent
theorems = {
    'MathCommons.lean': {'MathCommons.' + n for n in (
        'inEveryConvexSet', 'allPairings', 'orthogonal', 'determinantValue', 'witnessInHull')},
    'ForestCut.lean': {'MathCommons.ForestCut.' + n for n in (
        'telescoping', 'constantOnInterval', 'singleCutChange', 'totalChangeAtOpportunity',
        'sameEndpointsSameCrossing', 'cutHistoryUnique', 'crossingDeltaInjective', 'crossingDecisionUnique')},
}
allowed = {'propext', 'Classical.choice', 'Quot.sound'}
lean = sys.argv[1] if len(sys.argv) == 2 else 'lean'
version = (root / 'lean-toolchain').read_text().strip().split(':v')[1]
output = subprocess.check_output([lean, '--version'], cwd=root, text=True)
if not re.search(r'\bversion ' + re.escape(version) + r'[, )]', output):
    sys.exit(f'Expected Lean {version}, found {output.strip()}')
for source, expected in theorems.items():
    run = subprocess.run([lean, '-DwarningAsError=true', source], cwd=root, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print(run.stdout, end='')
    if run.returncode:
        sys.exit(run.returncode)
    seen = set()
    for name, deps in re.findall(r"'([\w.]+)' (?:does not depend on any axioms|depends on axioms:\s*\[([^\]]*)\])", run.stdout):
        used = {entry.strip() for entry in deps.split(',') if entry.strip()}
        if used - allowed:
            sys.exit(f'{name} depends on unapproved axioms: {sorted(used - allowed)}')
        seen.add(name)
    if expected - seen:
        sys.exit(f'Missing axiom reports: {sorted(expected - seen)}')
