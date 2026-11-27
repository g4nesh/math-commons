# Cold replay of P(154) · agent `claude-opus-5-5-fourier-154-replay`

**Verdict.** Every principal minor of the complex Fourier matrix of order 154 is nonzero, *conditional on* the correctness of Theorem 1.4 of Caragea–Lee–Malikiosis–Pfander (CLMP, [arXiv 2505.24326v3](https://arxiv.org/abs/2505.24326v3)). The theorem's finite hypothesis for (p, N') = (11, 14) is verified here by two independent exact computations. Evidence level: an exact finite computation plus a cited theorem; the theorem's proof was not re-derived.

**Independence.** This work was produced by Sonnet subagents coordinated by this agent for GitHub owner Raghavsk24, who is not the seed owner (g4nesh). Before this evidence was frozen, nobody here opened `order-154/seed.md`, `verify_seed.py`, `synthesis.md`, `prior-art.md`, family-14p material, or anyone's order-154 code. The two checkers were written without seeing each other.

## Argument

1. **Theorem 1.4** (`thm:inductionstep`). Let N = pN' be square-free with p prime, and suppose all principal minors of F_{N'} are nonzero in characteristic p. Then every N'-principal minor of F_N is nonzero, and hence so is every principal minor. [source-check.md](source-check.md) (a) quotes this from the LaTeX source. The theorem has no further hypothesis on p or N.
2. **Coverage.** A principal submatrix (K = L) is d-principal for every d (`defn:minors`), so every principal minor of F_154 falls under the conclusion. [source-check.md](source-check.md) (b).
3. **Choice.** 154 = 11 · 14 with p = 11 and N' = 14. Since ord₁₄(11) = 3, the residue field is F₁₃₃₁. Q(ω₁₄) has degree 6, so 11 splits into 2 primes. They correspond to ζ and ζ⁻¹, because −1 ∉ ⟨11⟩ = {1, 9, 11} in (Z/14)*.
4. **Meaning of "nonzero in characteristic p".** The proof needs p ∤ Norm(det Z_k), where Z_k are plain principal minors of F_{N'} evaluated at ζ = ω^p. So each minor must be nonzero modulo *every* prime above 11. Both ζ and ζ⁻¹ are checked. The paper's one-prime equivalence passes through Galois conjugates that need not be principal, so it is not relied on here. [source-check.md](source-check.md) (c, d).
5. **Finite check (Hyp(14, 11)).** For every S ⊆ Z₁₄ (all 2¹⁴ subsets, no quotient) and each root in {ζ, ζ⁻¹}, det(ζ^(ij))_{i,j∈S} ≠ 0 in F₁₃₃₁, using unnormalized entries and with the empty determinant equal to 1. That is 32 768 determinants, which exhausts the hypothesis as stated.
6. **Conclusion.** (5) satisfies the hypothesis of (1), and (2) gives every principal minor of F_154 ≠ 0 over ℂ. Nonvanishing is Galois-invariant, so the choice of primitive 154th root does not matter.

A failed finite check would only have refuted this lifting route, not the complex claim. See [reduction.md](reduction.md) §4.

## Computation

| | Checker A | Checker B |
|---|---|---|
| File | [checker_a/check_f1331.py](checker_a/check_f1331.py) | [checker_b/check_zech.py](checker_b/check_zech.py) |
| Field | F₁₁[x]/(x³ + x + 4), polynomial basis | F₁₁[x]/(x³ + 10x² + 10x + 5), Zech logarithms |
| Determinant | Gaussian elimination with inverses | Bareiss fraction-free, with cofactor cross-check |
| Main run | All 16 384 subsets × 2 roots | All 16 384 subsets × 2 roots, plus the quotient over 245 affine orbits |
| Result | 32 768 / 32 768 nonzero | 32 768 / 32 768 nonzero; quotient and full run agree |
| Transcript | [checker_a/run_transcript.txt](checker_a/run_transcript.txt) | [checker_b/run_transcript.txt](checker_b/run_transcript.txt) |

Counts by size are in each `counts.json`. Both roots have every subset of size k nonzero, for k = 0…14 (1, 14, 91, …, 14, 1).

Self-tests that pass in both checkers:
- field axioms and Frobenius (y¹³³¹ = y);
- the empty set gives 1;
- singletons give ζ^(i²) ≠ 0;
- M(ζ)M(ζ⁻¹) = 14I;
- zero detection on singular matrices, for example the order-7 root ζ² with S = {0, 7}, and 40 / 1024 zero minors at N = 10 over F₁₁ (checker B).

Checker A also tests the Jacobi complement identity (600 samples) and the translation scaling identity (300 samples). Checker B shows that all 84 affine maps preserve vanishing on all 32 768 cases.

**Environment.** Windows 11, Python 3.11.9, standard library only. Runtimes are about 154 s for checker A and 12 s for checker B.

**Commands.** Run from each checker folder:
- `python check_f1331.py`
- `python check_zech.py`

## Corrections and limits

- `reduction.md` says that singletons and size-1 minors equal 1. In fact they equal ζ^(i²), which is nonzero, and equals 1 only for i = 0. Both checkers test the correct statement. This does not affect the verdict.
- `reduction.md` was drafted from summarized paper text. Its open points 1–3 are resolved by [source-check.md](source-check.md), which reads the raw LaTeX: arXiv e-print sha256 `e70e507f60e2b222c83e9680c541fb8a0f9963071d29a04e2e3a532cf8648310`, downloaded 2026-11-27. Theorem and equation numbers there are LaTeX labels; the numbering in the compiled PDF was not checked.
- Not done:
  - the proof of Theorem 1.4 and of Lemma `mnprincipalpermutation` was not re-derived;
  - no direct complex cross-check at small N (14, 22, 42);
  - the CLMP Thm 4.6 smoke test was not run.
- The paper does not tabulate N' = 14. It lists 70 and 143 as its smallest unresolved square-free cases. Whether 154 is new is a question for `fourier-154-prior-art`.
