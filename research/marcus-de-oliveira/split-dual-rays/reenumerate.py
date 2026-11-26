"""Optional exact full enumeration with a locally built lrslib mini-GMP executable."""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import subprocess
import tempfile
import time

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('split_checker', HERE / 'check.py')
c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--lrs', required=True, help='path to a trusted lrsmgmp binary')
    args = parser.parse_args()
    binary = Path(args.lrs).resolve(strict=True)
    c.check()
    started = time.monotonic()
    run = subprocess.run([str(binary), str(HERE / 'bounded.ine')],
                         capture_output=True, text=True, timeout=60, check=True)
    c.a.require(not run.stderr.strip(), 'no enumerator stderr')
    with tempfile.TemporaryDirectory() as temp:
        output = Path(temp) / 'enumeration.txt'
        output.write_text(run.stdout)
        rays = c.validate_output(output)
    expected = {tuple(r) for r in json.loads((HERE / 'rays.json').read_text())}
    c.a.require(rays == expected, 'fresh complete enumeration agrees with every archived ray')
    print(json.dumps({'process_exit': run.returncode, 'rays': len(rays),
                      'exact_archived_set_equal': True,
                      'binary_sha256': hashlib.sha256(binary.read_bytes()).hexdigest(),
                      'elapsed_seconds': round(time.monotonic() - started, 3),
                      'tool': run.stdout.splitlines()[0]}, indent=2))


if __name__ == '__main__':
    main()
