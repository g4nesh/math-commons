// Canned role outputs for the Euler-polynomial fixture campaign. They stand in
// for model runs so the whole pipeline (gates, ingestion, certificates, Lean,
// publication) can be exercised without a model. The mathematics is classical
// (Euler, 1772); the campaign is marked as a fixture and never counted.

export const STATEMENT = `For every natural number $n \\ge 0$, the number $n^2 + n + 41$ is prime.`;

const statementLean = `namespace EulerPoly

/-- \`p\` is prime: at least 2, and its only divisors are 1 and itself. -/
def IsPrime (p : Nat) : Prop := 2 ≤ p ∧ ∀ d : Nat, d ∣ p → d = 1 ∨ d = p

/-- The source statement: \`n² + n + 41\` is prime for every natural number \`n\`. -/
def Claim : Prop := ∀ n : Nat, IsPrime (n ^ 2 + n + 41)

end EulerPoly
`;

const mainLean = `import EulerPoly.Statement

namespace EulerPoly

/-- C1: the value at n = 40 is 41². -/
theorem value_at_forty : 40 ^ 2 + 40 + 41 = 41 * 41 := by decide

/-- C3: the statement fails at n = 40. -/
theorem main : ¬ Claim := by
  intro h
  have h40 : 2 ≤ 40 ^ 2 + 40 + 41 ∧ ∀ d : Nat, d ∣ 40 ^ 2 + 40 + 41 → d = 1 ∨ d = 40 ^ 2 + 40 + 41 := h 40
  exact (h40.2 41 ⟨41, value_at_forty⟩).elim (fun h1 => absurd h1 (by decide)) (fun h2 => absurd h2 (by decide))

end EulerPoly
`;

const verifier = `"""Independent check of certificates/euler-forty: exact integer arithmetic only."""
import json
import sys


def check(data):
    n, value, factor = int(data["n"]), int(data["value"]), int(data["factor"])
    results = []
    results.append(("value equals n^2 + n + 41", value == n * n + n + 41))
    results.append(("factor divides value", value % factor == 0))
    results.append(("factor is nontrivial", 1 < factor < value))
    return results


def main():
    with open("certificate.json") as handle:
        data = json.load(handle)
    results = check(data)
    for label, ok in results:
        print(("PASS " if ok else "FAIL ") + label)
    corrupted = dict(data, factor=str(int(data["factor"]) + 1))
    assert not all(ok for _, ok in check(corrupted)), "negative control was accepted"
    print("PASS negative control rejected")
    if all(ok for _, ok in results):
        print("VERIFIED")
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())
`;

const proof = revision => `# Euler's polynomial is not prime-valued

<!-- digest -->
Complete disproof: at n = 40, n² + n + 41 = 1681 = 41², which is not prime.
Method: exhaustive search of n ≤ 1000 found the first failure; the value is certified
by an independent exact verifier and the argument is two lines. Nothing about the
density of prime values is claimed.
<!-- /digest -->

**Result kind:** \`complete\`
**Main claim:** \`C3\`

## Statement

The statement "for every natural number $n$, $n^2 + n + 41$ is prime" is false; it fails at $n = 40$.

## Claims

| ID | Statement | Status | Depends on | Evidence |
|---|---|---|---|---|
| C1 | $40^2 + 40 + 41 = 1681 = 41 \\cdot 41$ | computed | — | certificates/euler-forty |
| C2 | $1681$ is not prime | proved | C1 | divisor 41 with $1 < 41 < 1681$ |
| C3 | The universal statement fails | proved | C2 | instance $n = 40$ |

## Proofs

### C1. The value at forty

$40^2 + 40 + 41 = 1600 + 81 = 1681 = 41^2$; the certificate checks this in exact integer arithmetic.

### C2. Not prime

By C1, $41$ divides $1681$ and $1 < 41 < 1681$, so $1681$ has a divisor other than $1$ and itself.

### C3. The statement fails

The statement asserts primality for every $n$, in particular for $n = 40$; C2 contradicts this.

## Scope and limitations

Only the universal statement is refuted. Nothing is claimed about how many values are prime
(they are prime for $0 \\le n \\le 39$), nor about other prime-generating polynomials.

## Reproduction

\`python3 -I -B verify.py\` in \`certificates/euler-forty/\`.
${revision ? '\n## Revision notes\n\n- Addressed review-001 D1: stated the C2 divisor bounds explicitly.\n' : ''}`;

export function roleOutput(role, { runId, direction, attempt, certificate, mode, reviewsSoFar = 0 }) {
  switch (role) {
    case 'screen':
      return {
        files: {
          'SCREEN.md': `**Decision:** fixture (classical; used to exercise the harness)\n\n## Exact target\n\n${STATEMENT}\n\n## Prior-resolution search\n\n| Query | Where | Result |\n|---|---|---|\n| "n^2+n+41" prime | Google Scholar | Euler 1772; fails at n = 40 |\n| Euler prime-generating polynomial | Wikipedia | fails at n = 40 |\n| OEIS A005846 | OEIS | primes of the form n^2+n+41 |\n\n## Sources read\n\n- https://oeis.org/A005846\n- https://mathworld.wolfram.com/Prime-GeneratingPolynomial.html\n- https://en.wikipedia.org/wiki/Formula_for_primes\n\n## Search failures\n\nNone.\n\n## Significance\n\nNone: a fixture.\n\n## Full vs partial\n\nA single counterexample is a complete disproof.\n\n## Attacks\n\n- Exhaustive search over small n.\n`
        },
        answer: {
          verdict: mode === 'reject' ? 'prior-resolution' : 'fixture', summary: 'Known false since Euler (n = 40); kept as a harness fixture.',
          exactStatement: STATEMENT, source: { citation: 'L. Euler (1772), letter to Bernoulli', url: 'https://oeis.org/A005846', statementNumber: 'n/a' },
          priorResolutionSearch: {
            queries: [{ query: 'n^2+n+41 prime', where: 'Google Scholar', result: 'fails at n = 40' }, { query: 'Euler prime-generating polynomial', where: 'Wikipedia', result: 'fails at n = 40' }, { query: 'A005846', where: 'OEIS', result: 'sequence of prime values' }],
            sources: [{ url: 'https://oeis.org/A005846', title: 'OEIS A005846', finding: 'prime values listed' }, { url: 'https://mathworld.wolfram.com/Prime-GeneratingPolynomial.html', title: 'MathWorld', finding: 'composite at n = 40' }, { url: 'https://en.wikipedia.org/wiki/Formula_for_primes', title: 'Formula for primes', finding: 'composite at n = 40' }],
            failures: []
          },
          significance: 'None. This is a classical fact used only to test every stage of the harness end to end without a model.',
          fullVsPartial: 'One counterexample refutes the universal statement completely.',
          attacks: [{ name: 'search', idea: 'evaluate small n', firstTest: 'n ≤ 1000' }],
          priorResolution: mode === 'reject' ? 'Euler 1772: composite at n = 40' : null
        }
      };
    case 'map':
      return {
        files: { 'MAP.md': '# Map\n\n## Objects and conventions\n\n$n \\in \\mathbb{N}$ including 0; prime means at least 2 with no nontrivial divisor.\n\n## Smallest cases\n\nValues 41, 43, 47, … are prime for n ≤ 39.\n\n## Known tools\n\nFactor $n^2+n+41$ at $n \\equiv 0, 40 \\pmod{41}$.\n\n## Lenses\n\nAdversarial computation; modular structure.\n\n## Where it should break\n\nAt n = 40 and n = 41, where 41 divides the value.\n' },
        answer: {
          summary: 'Two directions: exhaustive search and modular structure.',
          directions: [
            { id: 'D1', title: 'Exhaustive search', lens: 'adversarial computation', kind: 'adversarial', hypothesis: 'a small counterexample exists', falsificationTest: 'all n ≤ 1000', milestone: 'an exact counterexample', budgetTokens: 200000 },
            { id: 'D2', title: 'Divisibility by 41', lens: 'modular structure', kind: 'constructive', hypothesis: '41 divides the value when n ≡ 0 or 40 mod 41', falsificationTest: 'check n = 40, 41', milestone: 'a family of composite values', budgetTokens: 200000 }
          ]
        }
      };
    case 'adversary':
      return {
        files: {
          [`attempts/${direction}/${attempt}.md`]: `**Outcome:** counterexample · **Status:** promising\n\n<!-- digest -->\nFirst failure at n = 40: 1681 = 41². Exhaustive for n ≤ 1000.\n<!-- /digest -->\n\n## Search\n\nAll n in [0, 1000]; first composite value at n = 40.\n`,
          'certificates/euler-forty/certificate.json': JSON.stringify({ n: '40', value: '1681', factor: '41' }, null, 2) + '\n',
          'certificates/euler-forty/search.py': 'for n in range(1001):\n    v = n * n + n + 41\n    if any(v % d == 0 for d in range(2, int(v ** 0.5) + 1)):\n        print(n, v)\n        break\n',
          'certificates/euler-forty/manifest.json': JSON.stringify({ name: 'euler-forty', claim: 'C1: 40^2 + 40 + 41 = 1681 has the nontrivial divisor 41', description: 'Counterexample at n = 40', notClaimed: 'Nothing about other n', generator: { path: 'search.py', command: ['python3', 'search.py'] }, verifier: { path: null, command: null } }, null, 2) + '\n'
        },
        answer: { summary: 'n = 40 gives 1681 = 41².', file: `attempts/${direction}/${attempt}.md`, certificates: ['euler-forty'], attempt: { claim: 'The statement fails at n = 40', method: 'exhaustive search n ≤ 1000', outcome: 'counterexample', status: 'promising', result: '1681 = 41 · 41', failurePoint: null, evidence: ['certificates/euler-forty'], next: 'verify the certificate independently' } }
      };
    case 'explore':
      return {
        files: { [`attempts/${direction}/${attempt}.md`]: `**Outcome:** partial · **Status:** open\n\n<!-- digest -->\n41 divides the value whenever n ≡ 0 or 40 (mod 41).\n<!-- /digest -->\n\nAt n = 41k, the value is 41(41k² + k + 1).\n` },
        answer: { summary: 'A family of composite values.', file: `attempts/${direction}/${attempt}.md`, certificates: [], attempt: { claim: '41 | n²+n+41 for n ≡ 0, 40 mod 41', method: 'modular arithmetic', outcome: 'partial', status: 'open', result: 'infinitely many composite values', failurePoint: null, evidence: [], next: 'none needed once D1 is drafted' } }
      };
    case 'verify-certificate':
      return {
        files: {
          [`certificates/${certificate}/verify.py`]: verifier,
          [`certificates/${certificate}/VERIFIER.md`]: '# Verifier checklist\n\n1. value = n² + n + 41\n2. factor divides value\n3. 1 < factor < value\n',
          [`certificates/${certificate}/manifest.json`]: JSON.stringify({ name: 'euler-forty', claim: 'C1: 40^2 + 40 + 41 = 1681 has the nontrivial divisor 41', description: 'Counterexample at n = 40', notClaimed: 'Nothing about other n', generator: { path: 'search.py', command: ['python3', 'search.py'] }, verifier: { path: 'verify.py', command: ['python3', '-I', '-B', 'verify.py'] } }, null, 2) + '\n'
        },
        answer: { summary: 'Verifier written; passes.', verdict: 'pass', command: ['python3', '-I', '-B', 'verify.py'], checks: ['value', 'divisibility', 'nontrivial factor'] }
      };
    case 'draft':
      return { files: { 'PROOF.md': proof(false) }, answer: { summary: 'Complete disproof at n = 40.', kind: 'complete', main: 'C3', polarity: 'disprove', claims: [{ id: 'C1', status: 'computed' }, { id: 'C2', status: 'proved' }, { id: 'C3', status: 'proved' }], files: ['PROOF.md'] } };
    case 'review': {
      const changes = mode === 'review-changes-once' && reviewsSoFar === 0;
      return {
        files: { [`reviews/work-${runId}/REPORT.md`]: `**Verdict:** ${changes ? 'changes' : 'pass'}\n\n## Target match\n\nExact.\n\n## Claims checked\n\nC1, C2, C3.\n\n## Defects\n\n${changes ? 'D1 (blocking): C2 does not state 1 < 41 < 1681.' : 'None.'}\n\n## Independent checks\n\nRe-ran verify.py: VERIFIED.\n\n## Status boundary\n\nNone beyond the certificate re-run.\n` },
        answer: {
          verdict: changes ? 'changes' : 'pass', summary: changes ? 'C2 needs explicit divisor bounds.' : 'Every claim checked; the counterexample is exact.',
          checkedClaims: ['C1', 'C2', 'C3'],
          defects: changes ? [{ id: 'D1', severity: 'blocking', location: 'C2', description: 'C2 must state 1 < 41 < 1681.', suggestedRepair: 'State the bounds.' }] : [],
          independentChecks: [{ what: 'certificate', command: 'python3 -I -B verify.py', outcome: 'VERIFIED' }],
          statusBoundary: 'Certificate checked by re-running its verifier.', kindAccurate: true
        }
      };
    }
    case 'repair':
      return { files: { 'PROOF.md': proof(true) }, answer: { summary: 'D1 addressed.', addressed: [{ defect: 'review-001#D1', change: 'bounds stated' }], unresolved: [], kind: 'complete' } };
    case 'formalize-statement':
      return {
        files: {
          'formal/lakefile.toml': 'name = "eulerpoly"\nversion = "0.1.0"\ndefaultTargets = ["EulerPoly"]\n\n[[lean_lib]]\nname = "EulerPoly"\n',
          'formal/lean-toolchain': 'leanprover/lean4:v4.24.0\n',
          'formal/EulerPoly.lean': 'import EulerPoly.Main\n',
          'formal/EulerPoly/Statement.lean': statementLean,
          'formal/EulerPoly/Main.lean': 'import EulerPoly.Statement\n',
          'STATEMENT_REVIEW.md': '# Statement review\n\n## Formalization notes\n\n| Source condition | Lean representation | Note |\n|---|---|---|\n| every natural number n ≥ 0 | `∀ n : Nat` | Lean ℕ includes 0 |\n| prime | `IsPrime p := 2 ≤ p ∧ ∀ d, d ∣ p → d = 1 ∨ d = p` | standard |\n| n² + n + 41 | `n ^ 2 + n + 41` | literal |\n'
        },
        answer: { summary: 'Statement formalized.', status: 'done', files: ['formal/EulerPoly/Statement.lean'], lean: { modules: ['EulerPoly.Main'], claim: 'EulerPoly.Claim', main: 'EulerPoly.main', namespace: 'EulerPoly', endpoints: ['EulerPoly.value_at_forty'] }, gaps: [] }
      };
    case 'audit-statement':
      return {
        append: { 'STATEMENT_REVIEW.md': '\n## Audit\n\n**Verdict:** faithful\n\n| Source condition or quantity | Literal Lean representation | Match | Note |\n|---|---|---|---|\n| for every natural number n | `∀ n : Nat` | exact | includes 0 as the source does |\n| prime | `IsPrime` | exact | 0 and 1 excluded by `2 ≤ p` |\n| n² + n + 41 | `n ^ 2 + n + 41` | exact | |\n' },
        answer: { verdict: 'faithful', summary: 'Every condition matches literally.', mapping: [{ source: 'for every n', lean: '∀ n : Nat', match: 'exact', note: '' }], defects: [] }
      };
    case 'formalize':
      return { files: { 'formal/EulerPoly/Main.lean': mainLean }, answer: { summary: 'Proved.', status: 'done', files: ['formal/EulerPoly/Main.lean'], lean: { modules: ['EulerPoly.Main'], claim: 'EulerPoly.Claim', main: 'EulerPoly.main', namespace: 'EulerPoly', endpoints: ['EulerPoly.value_at_forty'] }, gaps: [] } };
    default:
      throw new Error(`fixture has no output for role ${role}`);
  }
}
