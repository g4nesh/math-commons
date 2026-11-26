Evidence level: PROOF-READ (Theorem 1, line by line, standard potential theory supplied where the text is terse); EXACT (the algebraic facts of the "height exception" and application paragraphs, replayed with own code); NUMERIC (one illustration of the pole step).

# Audit: final_referee_bounded_bad_harmonic_lemma_20261023

PDF: `proofs/manin_transfer_20261021/final_referee_bounded_bad_harmonic_lemma_20261023.pdf` (2 pages). Batch F4-common-lemmas-b.
Replay: `replay/final_referee_bounded_bad_harmonic_lemma_20261023__check.py` (output `..._check_output.txt`, 4/4 EXACT pass, 1 NUMERIC pass).

## Deviations and limits

- Theorem 1 is an analytic lemma with no machine-checkable content; it was read line by line. Standard facts used without re-proof: a continuous function satisfying the local sub-mean-value inequality is subharmonic; a subharmonic function bounded above near an isolated point extends subharmonically across it; a subharmonic function on a compact connected Riemann surface is constant. These are textbook facts (e.g. Ransford, Potential Theory in the Complex Plane, Thms 2.4.1, 3.6.1 and the maximum principle), recorded as PROOF-READ, not as dependencies.
- The document is a common lemma; whether each application verifies hypotheses 1 to 3 (in particular boundedness of d near every puncture) is the job of the auditors of the citing documents. This report lists what each application must check.
- No supplement code or data was used.

## A. Statements

**Theorem 1 (p. 1).** Let X be a compact connected Riemann surface, P ⊂ X finite, H_1..H_N meromorphic on X, holomorphic on X \ P. Fix real η ≠ 0, g_i = −η Im H_i (harmonic on X \ P), and L = Σ λ_i g_i (real λ_i, no constant term). Let D be a connected component of the open set {x ∈ X \ P : g_i(x) > 0 ∀i, L(x) < 0}. Assume
1. at every point x ∈ X \ P on the boundary of D: if g_i(x) ≥ 0 ∀i and some g_i(x) = 0, then L(x) ≥ 0;
2. there are d : D → [0, ∞) and M ≥ 0 with −L ≤ M d on D, and d is bounded on D ∩ U_p for a punctured neighbourhood U_p of every p ∈ P;
3. not all H_i are constant.

Then D = ∅. ("Ordinary boundary point" in the text means a boundary point of D in X \ P; the strictly supported negative set is {g > 0, L < 0}.)

**Application remarks (p. 2), not theorems.** (a) If c = Σ_k λ_k f_k + d f_bad with minor f_k, λ_k ≥ 0 and d ≥ 0 on D, then L = Φ_c(O) ≥ d Φ_{f_bad}(O), so M = max(0, −Φ_{f_bad}(O)) works; (b) if d ≤ d + g_j with d + g_j bounded near p, hypothesis 2 holds at p; (c) "height exception": if all gaps are constant on the curve but x = f κ + v_0 with f nonconstant, then ∇C(v_0) = 0 and H_κ v_0 = 0, so v_0 lies in the span of single-entry vectors and the cost is Birkhoff-safe; a constant full matching vector must be excluded separately.

## B. Dependencies

| dependency | kind | status |
|---|---|---|
| removable singularities / maximum principle for subharmonic functions | standard analysis | PROOF-READ (textbook) |
| ker H_κ = span{n_{i,j}}, ∇C(κ) = 0, X ⊂ Sing C (height exception) | paper (3.10)-(3.12), Appendix A | EXACT (own replay, this script and moc4lib) |
| Birkhoff-von Neumann theorem (height exception) | external classical | CITED |
| Φ_{n_{I,J}}(O) = det O[I,J]^2 ≥ 0 (application (a)) | paper (4.2) | EXACT (moc4lib, plus this script on 10 exact O) |

No contact case, branch document, Appendix B or Corollary B.3 is invoked. The "positive-gap compensation" remark (b) is the trivial inequality d ≤ d + g_j on D, not Corollary B.3.

## C. Steps

| # | claim | kind | status | how checked |
|---|---|---|---|---|
| 1 | v := −L on D, 0 elsewhere; at a boundary point q ∈ X \ P of D, L(q) = 0 | analytic | PROOF-READ | By continuity g_i(q) ≥ 0, L(q) ≤ 0. If all g_i(q) > 0 and L(q) < 0, q lies in the open set and its component is D, so q ∈ D, a contradiction. Hence some g_i(q) = 0 (then hyp. 1 gives L(q) ≥ 0) or L(q) = 0. Either way L(q) = 0. |
| 2 | v is subharmonic on X \ P | analytic | PROOF-READ (argument supplied by auditor) | v is continuous by step 1, harmonic on D, 0 near points outside cl D, and at points with v = 0 the sub-mean-value inequality holds because v ≥ 0. The text's "pasting lemma" is this. |
| 3 | v extends subharmonically over P | analytic | PROOF-READ | v ≤ M sup d near p (hyp. 2), v = 0 off D; finite sets are polar, so a function bounded above near p extends. |
| 4 | v constant > 0 on X, hence D = X \ P | analytic | PROOF-READ | Maximum principle on compact connected X; v > 0 somewhere since D ≠ ∅; then v > 0 everywhere, and v > 0 only on D. |
| 5 | every g_i > 0 on X \ P forces H_i holomorphic at P | analytic | PROOF-READ; NUMERIC illustration | At a pole of order k, H = c z^{-k}(1 + O(z)) and Im H takes both signs on small circles. So g_i takes both signs on a punctured disc. |
| 6 | all H_i constant, contradicting hyp. 3 | analytic | PROOF-READ | Holomorphic functions on a compact connected surface are constant. |
| 7 | application (a): Φ_c(O) ≥ d Φ_{f_bad}(O) | inequality | EXACT on instances; PROOF-READ in general | Φ_{f_k}(O) = det^2 ≥ 0 (paper (4.2)). The script checks it for all 48 non-minor f_bad on 10 exact rational O. Nonnegativity of the λ_k on D is NOT part of the lemma; each application must supply it. |
| 8 | application (b): d ≤ d + g_j bounded suffices | inequality | PROOF-READ | g_j > 0 on D. |
| 9 | height exception: ∇C(v_0 + fκ) = ∇C(v_0) + f H_κ v_0, ker H_κ = span{n_ij} | identity | EXACT | Symbolic in all 14 coordinates; ker H_κ has dim 10 and equals the span of the 16 single-entry vectors. |
| 10 | constant supported single-entry cost is Birkhoff-safe | inequality | PROOF-READ (+ Birkhoff CITED) | Φ = Σ a_ij O_ij^2 and (O_ij^2) is doubly stochastic (EXACT on instances). |

## D. Falsification

Theorem 1 is a purely analytic statement, and no counterexample search makes sense for it. Its weakest points are the hypotheses an application has to supply. I tried the two places where the hypotheses could silently fail:
- If a gap is identically zero on the curve, the set {g_i > 0 ∀i} is empty and the lemma says nothing. Contact gaps must therefore be excluded from the H_i, as the application paragraph does ("noncontact ... coordinate functions").
- Hypothesis 2 is needed at EVERY p ∈ P, including the punctures where a spectral parameter goes to infinity and those where two parameters collide. The lemma gives no help there. Each citing document must list all punctures.
NUMERIC: on 200 random Laurent heads with pole order 1 to 4 (seed 1), Im H took both signs on |z| = 1e-3.

## E. Verdict

**confirmed.** Theorem 1 is correct as stated. Every step is PROOF-READ, and the algebraic remarks are EXACT. Birkhoff's theorem is used only in the height-exception remark, not in Theorem 1.

Presentation notes (not gaps): "ordinary boundary point" and "strictly supported negative set" are not defined, but the meaning is clear from the proof. The sentence beginning "Equivalently, D is a union of components ..." is not needed. Obligations on every application, to be checked by the auditors of the citing documents: (i) the H_i are the noncontact gaps, meromorphic on the compactified curve and holomorphic off the listed P; (ii) a decomposition with λ_k ≥ 0 and d ≥ 0 on D; (iii) d bounded near every puncture, all punctures listed; (iv) nonconstancy (hyp. 3), or the height exception with the constant-vector case excluded separately.
