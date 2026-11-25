# Exact finite-witness certificates

This directory verifies **one specified instance** of the original Marcus–de Oliveira convex-hull statement using Gaussian-rational arithmetic. It does not prove the universal conjecture, validate an arbitrary floating-point search, or supply independent peer review. The tests and reproduced results here were produced by agents under the same GitHub owner; another owner's review is pending.

For input spectra `a,b` and a unitary `U`, the program constructs

$$
A=\operatorname{diag}(a),\quad B=U\operatorname{diag}(b)U^*,\quad d=\det(A+B),
\qquad z_\sigma=\prod_i(a_i+b_{\sigma(i)}).
$$

It checks both `U U* = I` and `U* U = I` exactly, constructs `B`, verifies both matrices' normality, and computes `d` by exact Gaussian elimination. It enumerates **all** `n!` permutation points. Membership certificates contain at most three nonnegative rational weights summing to one. Outside certificates contain a rational real linear functional whose value at `d` exceeds its value at every permutation point by a strictly positive rational gap.

## Run

Python 3.9+; standard library only. From the repository root:

```sh
python3 -m unittest discover -s research/marcus-de-oliveira/certificates -p 'test_*.py' -v
python3 research/marcus-de-oliveira/certificates/verify.py research/marcus-de-oliveira/certificates/fixtures/goyal-2026.json
python3 research/marcus-de-oliveira/certificates/verify.py research/marcus-de-oliveira/certificates/fixtures/zhang-ma-2026.json
python3 research/marcus-de-oliveira/certificates/verify.py research/marcus-de-oliveira/certificates/fixtures/complex-householder-n4.json
```

To save a new result, add `--output /path/to/result.json`. Exit status `0` means the exact computation succeeded; inspect `membership` for `inside` or `outside`. Exit status `2` means invalid input or an uncertified computation. An `outside` result would require independent reproduction and review before any claim of a counterexample.

## Witness format

See [witness.schema.json](witness.schema.json). Each complex number is `[real, imaginary]`, with both components **strings**. Use integer or `p/q` strings, with a positive denominator. The two spectra have the same length `1 ≤ n ≤ 6`, and `U` is an `n × n` matrix. Each integer component has at most 64 decimal digits; files are limited to 1 MiB. Duplicate JSON keys, unknown fields, JSON numbers in complex entries, decimals, exponents, zero/negative denominators, mismatched shapes, and nonunitary matrices are rejected.

```json
{
  "schemaVersion": 1,
  "name": "singular one-dimensional example",
  "a": [["1", "1"]],
  "b": [["-1", "-1"]],
  "U": [[["0", "1"]]],
  "expectedDeterminant": ["0", "0"]
}
```

Optional `name`, `source`, and `notes` are strings. Optional `expectedDeterminant` is an assertion checked **after recomputation**; it never replaces the computed determinant. `source` is metadata and is never fetched or executed.

All output complex pairs use the same exact string representation. Permutations use **one-based indices**: `[2,1,3]` means `b₂,b₁,b₃` are paired with `a₁,a₂,a₃`. Duplicate permutation points remain in `permutationPoints`; hull construction deduplicates coordinates and uses one representative permutation for each certificate weight.

## Reproduced instances

The first two fixtures transcribe spectra and matrices from Section 4 of [Kovacec, DMUC 26-11, v2, November 4, 2026](https://www.mat.uc.pt/preprints/ps/p2611.pdf#page=8). Those examples concern failure of a stronger polygon-enclosure statement. The computations below independently establish membership in the **original convex hull** for the supplied inputs; this checker does not construct that stronger polygon.

| Input | Recomputed determinant | Original convex-hull membership |
| --- | --- | --- |
| [Goyal input](fixtures/goyal-2026.json) | `2/3 + (2/9)i` | Inside; [exact result](results/goyal-2026.json) |
| [Zhang–Ma input](fixtures/zhang-ma-2026.json) | `(2864 + 1840i)/49` | Inside; [exact result](results/zhang-ma-2026.json) |
| [Complex Householder, n=4](fixtures/complex-householder-n4.json) | `-880/27 + (800/81)i` | Inside; [exact result](results/complex-householder-n4.json) |

Two compact certificates from this implementation are

$$
\frac23+\frac29i
=\frac{13}{27}(-i)+\frac13(2+i)+\frac5{27}(2i),
$$

and

$$
\frac{2864+1840i}{49}
=\frac4{49}(-130+10i)+\frac{36}{49}(24-20i)+\frac9{49}(280+280i).
$$

Each point on the right is one of the fully enumerated permutation points in its result file. The n=4 fixture is a constructed regression case: `v=(1,i,1+i,2−i)` has `v* v=9`, and `U=I−2vv*/9`. Its genuinely complex off-diagonal entries test conjugation, not only real orthogonal arithmetic. No novel general case is established by this finite example.

## Why the certificate check is sufficient for one input

The monotone-chain algorithm orders exact rational plane coordinates and removes non-extreme collinear points. A single extreme point is handled by equality. A segment is handled by exact collinearity and an interval parameter. For a two-dimensional convex hull, counterclockwise edge orientations decide membership. A fan triangulation then supplies barycentric weights for an inside point; boundary weights may be zero and are omitted.

A **separate direct certificate check**, which does not call hull construction, validates nonnegative weights, their sum, their support in the complete permutation-point set, and exact reconstruction of `d`. For an outside result it checks every permutation point against the proposed bound and verifies a strictly positive exact gap. Thus a hull-construction mistake cannot silently pass a mathematically invalid generated certificate through that direct check. This is an implementation safeguard, not a separately authored proof assistant or independent review.

The normality and spectrum hypotheses follow from the construction once exact unitarity is checked: a diagonal matrix is normal, and unitary conjugation preserves normality and eigenvalues. The additional normality multiplications are redundant structural regressions. Singularity of `A+B` is allowed; exact zero pivots return determinant zero rather than applying an invalid inverse.

## Validation and limitations

The 16 regression tests cover:

- Gaussian-rational arithmetic and exact conjugation; elimination against a separately implemented Leibniz determinant on 32 deterministic complex matrices of dimensions 1–4.
- Both published n=3 instances, including their complete permutation-point sets, and the genuinely complex n=4 instance.
- Singular sums, repeated spectra, point hulls, segment hulls, vertices, boundary points, and strict outside synthetic points.
- A 169-point rational grid against an independently known square-membership condition; an outside gap as small as `1/10^30` in its coordinate.
- Rejection of corrupted convex weights/separation gaps and malformed, nonunitary, or unsupported witness inputs. The separately callable direct checker rejects floats and booleans in every scalar field; a regression covers a false singleton equality caused by rounding beyond `2^53`.
- Dimension six with all 720 permutations.

The **outside geometry tests are synthetic hull tests**, not MDO counterexamples. All three committed matrix fixtures certify inside. The tool supports `Q(i)` inputs only. It neither encloses arbitrary real/complex inputs with intervals nor proves the correctness of decimal approximations. In particular, the sibling numerical search's stored binary64 matrices cannot be converted into exact witnesses by merely spelling the rounded entries as fractions: exact unitarity will generally fail. Use an exact rational unitary construction or a separate rigorous interval argument.

The dimension and input limits bound enumeration and input size; they are not a universal theorem, a formal runtime guarantee, or coverage of every rational witness. Saved JSON files are deterministic replay artifacts; tests compare current recomputation with them. Keep an immutable source commit and fixture hashes with any external certificate review.
