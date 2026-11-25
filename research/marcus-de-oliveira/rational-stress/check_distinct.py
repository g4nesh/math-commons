#!/usr/bin/env python3
"""Count distinct exact coordinate matrix pairs, independently of case labels.

This does not quotient by unitary similarity, permutations or MDO symmetries.
It hashes canonical A's diagonal and B=U diag(b) U*, excluding metadata and U.
"""
import hashlib
import json
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import screen as S


def main():
    start = time.process_time()
    seen, records, reports = {}, [], {}
    for relative in ('report.json','holdout/report.json'):
        path = HERE/relative
        reports[relative] = hashlib.sha256(path.read_bytes()).hexdigest()
        report = json.loads(path.read_text())
        for row in report['records']:
            witness = S.make_witness(row['spec'])
            if S.digest(witness) != row['witnessSha256']:
                raise ArithmeticError('Input regeneration mismatch')
            a,b = [list(map(S.V.parse_complex,witness[key])) for key in ('a','b')]
            U = [list(map(S.V.parse_complex,values)) for values in witness['U']]
            B = S.V.multiply(S.V.multiply(U,S.V.diagonal(b)),S.V.adjoint(U))
            mathematical_pair = {'a':[value.encode() for value in a],
                                 'B':[[value.encode() for value in values] for values in B]}
            key = S.digest(mathematical_pair)
            records.append({'id':row['spec']['id'],'report':relative,'matrixPairSha256':key,
                            'duplicatesEarlierId':seen.get(key)})
            seen.setdefault(key,row['spec']['id'])
    result = dict(schemaVersion=1,scope='Distinct exact ordered coordinate matrices (A,B), not symmetry-equivalence classes.',
                  cases=len(records),distinctMatrixPairs=len(seen),duplicateMatrixPairs=len(records)-len(seen),
                  reportSha256=reports,cpuSeconds=time.process_time()-start,records=records)
    S.write_json(HERE/'distinct-inputs.json',result)
    print(json.dumps({key:result[key] for key in ('cases','distinctMatrixPairs','duplicateMatrixPairs','cpuSeconds')}))


if __name__ == '__main__': main()
