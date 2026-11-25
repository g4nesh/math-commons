#!/usr/bin/env python3
"""Replay saved approximate samples; does not generate a new search."""
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from search import analyze_case


def decode(value):
    if not isinstance(value, list) or len(value) != 2 or not all(isinstance(part, (int, float)) for part in value):
        raise ValueError('Expected approximate [real number, imaginary number], not rational-string exact schema.')
    return complex(*value)


def replay(path):
    report = json.loads(path.read_text())
    if report.get('schema') != 'mdo.approximate-complex-pair-search.v1' or report.get('approximate') is not True:
        raise ValueError('This script accepts only the explicit approximate search schema.')
    mismatches = []
    for case in report['cases']:
        if 'error' in case:
            raise ValueError(f"Recorded generation error: {case['id']}")
        fresh = analyze_case(list(map(decode, case['a'])), list(map(decode, case['b'])), [[decode(value) for value in row] for row in case['U']])
        if fresh != case['analysis']:
            keys = [key for key in fresh if fresh[key] != case['analysis'].get(key)]
            mismatches.append({'id': case['id'], 'different_fields': keys})
    result = {'schema': 'mdo.saved-sample-replay.v1', 'input': str(path), 'cases_replayed': len(report['cases']), 'new_search_samples': 0, 'all_analysis_fields_identical': not mismatches, 'mismatches': mismatches, 'note': 'Bitwise numeric agreement is expected on the original Python/platform. Other libm/platform versions may differ and require review; a mismatch is not a mathematical counterexample.'}
    print(json.dumps(result, indent=2))
    return not mismatches


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report', type=Path, nargs='?', default=Path(__file__).with_name('report.json'))
    arguments = parser.parse_args()
    sys.exit(0 if replay(arguments.report) else 1)
