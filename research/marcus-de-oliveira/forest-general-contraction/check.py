"""Fresh Lean 4.24 replay, exact statement comparison and transitive axiom audit.

Run: python3 -I -B research/marcus-de-oliveira/forest-general-contraction/check.py /path/to/lean
Only temporary .olean files are written. No network, Mathlib, or cached builds.
"""

import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
FORMAL = ROOT / 'formal'
MODULE = 'ForestGeneralContraction'
ALLOWED = {'propext', 'Classical.choice', 'Quot.sound'}
REPORT = re.compile(r"'([\w.]+)' (?:does not depend on any axioms|depends on axioms:\s*\[([^\]]*)\])")


def sources():
    ordered, seen = [], set()

    def visit(name):
        if name in seen:
            return
        seen.add(name)
        path = FORMAL / (name + '.lean')
        if not path.is_file():
            raise RuntimeError(f'Unrecognized local dependency: {name}')
        for dependency in re.findall(r'^import\s+(\w+)\s*$', path.read_text(), re.M):
            if dependency != 'Std':
                visit(dependency)
        ordered.append(path)

    visit(MODULE)
    return ordered


def replay(lean):
    version = subprocess.check_output([lean, '--version'], text=True, cwd=FORMAL)
    if not re.search(r'\bversion 4\.24\.0[, )]', version):
        raise RuntimeError(f'Expected Lean 4.24.0, got {version.strip()}')
    manifest = json.loads((HERE / 'declarations.json').read_text())
    names = [manifest['namespace'] + '.' + name for name in manifest['declarations']]
    if len(names) != manifest['count'] or len(set(names)) != len(names):
        raise RuntimeError('Invalid declaration inventory')
    reports = {}
    paths = sources()
    hashes = {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}
    with TemporaryDirectory(prefix='math-commons-general-contraction-') as build:
        env = {**os.environ, 'LEAN_PATH': build}
        for path in paths:
            run = subprocess.run([lean, '-DwarningAsError=true', '-o', str(Path(build, path.stem + '.olean')), path.name],
                                 cwd=FORMAL, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            if run.returncode:
                raise RuntimeError(run.stdout)
            found = {}
            for name, dependencies in REPORT.findall(run.stdout):
                used = sorted(part.strip() for part in dependencies.split(',') if part.strip())
                if set(used) - ALLOWED:
                    raise RuntimeError(f'Unapproved axioms for {name}: {used}')
                if name in reports or name in found:
                    raise RuntimeError(f'Duplicate axiom report: {name}')
                found[name] = used
            expected = set(re.findall(r'^#print axioms (\S+)\s*$', path.read_text(), re.M))
            if len(found) != len(expected):
                raise RuntimeError(f'Missing axiom reports in {path.name}')
            reports.update(found)
        actual = {name: deps for name, deps in reports.items() if name.startswith(manifest['namespace'] + '.')}
        if set(actual) != set(names):
            raise RuntimeError('New declaration names do not match the manifest')
        statement_file = Path(build, 'Statements.lean')
        statement_file.write_text('import ' + MODULE + '\nset_option format.width 100\n' +
                                  '\n'.join('#check ' + name for name in names) + '\n')
        run = subprocess.run([lean, '-DwarningAsError=true', str(statement_file)], cwd=FORMAL, env=env,
                             text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if run.returncode:
            raise RuntimeError(run.stdout)
        statements = run.stdout
    return {'leanVersion': '4.24.0', 'fatalWarnings': True, 'moduleCount': len(paths),
            'transitiveDeclarationCount': len(reports), 'newDeclarationCount': len(actual),
            'sourceSha256': hashes, 'axiomReports': actual}, statements


def main():
    lean = sys.argv[1] if len(sys.argv) == 2 else 'lean'
    report, statements = replay(lean)
    if report != json.loads((HERE / 'audit.json').read_text()):
        raise RuntimeError('Source hashes or axiom audit differ from the published manifest')
    if statements != (HERE / 'statements.txt').read_text():
        raise RuntimeError('Exact Lean theorem statements differ from the published snapshot')
    print(f"PASS: {report['moduleCount']} freshly built modules; {report['transitiveDeclarationCount']} axiom reports")
    print(f"PASS: {report['newDeclarationCount']} exact public declarations and source hashes; no unapproved axioms")


if __name__ == '__main__':
    main()
