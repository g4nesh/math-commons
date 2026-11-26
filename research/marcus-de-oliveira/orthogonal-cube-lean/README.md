# Lean proof of the unit-cube quadratic lemma

[`formal/OrthogonalCube.lean`](../../../formal/OrthogonalCube.lean) proves Claude's **Lemma L** from [PR48, commit `41ca275e7354043979278d72f8db897b9e2905a0`](https://github.com/g4nesh/math-commons/tree/41ca275e7354043979278d72f8db897b9e2905a0/research/marcus-de-oliveira/star-orthogonal). The Lemma L argument is unchanged from `ec52b46c2e2a9ed21947fcf6d91a27b1c56d6658`; the newer source adds surrounding results, which are not formalized here.

For nonnegative parameters \(p,q,r,u\) with \(p+q=r+u\), define

\[
Q(a,b,c,d)=pr(a^2+c^2)+qu(b^2+d^2)-2pqab-2rucd.
\]

For every \((a,b,c,d)\in[0,1]^4\), the checked theorem `MathCommons.OrthogonalCube.lemmaL` states

\[
0\le Q(a,b,c,d)\quad\text{or}\quad0\le Q(1-a,1-b,1-c,1-d).
\]

`lemmaLMax` gives the equivalent nonnegative-maximum formulation. Its `orderedMax` is defined as `if x ≤ y then y else x` in the same lawful linear order. Neither final theorem assumes F1–F4, the product ordering \(pq\le ru\), positivity of any parameter, nor a strict interior point. All zero parameters and the entire closed cube are included.

## Exact domain and scope

The domain is an arbitrary type `R` carrying `Lean.Grind.Field R`, `LE R`, `LT R`, `Std.LawfulOrderLT R`, `Std.IsLinearOrder R`, and `Lean.Grind.OrderedRing R`. These are Lean's standard field laws and compatible lawful linear-order laws. The theorem is a general ordered-field statement; it is not a rational-only theorem, a floating-point check, or a theorem with auxiliary polynomial signs added to its hypotheses. It uses neither completeness nor square roots.

This file imports only `Std`. It does not import or instantiate Mathlib's concrete `Real` type. A real-field instance satisfying these standard laws is covered by the generic statement; connecting a particular external real-number implementation is separate library integration. The scoped checker also compiles a rational specialization as an inhabited-instance sanity check, not as a replacement for the generic proof.

There is **no formal matrix, SO(3), SO(4), Hodge, squared-minor, or PMSH reduction in this module**. Those bridges must be checked separately before using this lemma to claim a matrix result. Claude's newer surrounding Theorem W is not included. This is same-owner work; a different GitHub owner must perform the ledger review.

## Proof organization

Let \(q_1(x,y)=prx^2-2pqxy+quy^2\) and \(q_2(x,y)=prx^2-2ruxy+quy^2\). In the case \(pq\le ru\), the proof derives the following facts rather than assuming them.

| Step | Checked declarations | Role |
| --- | --- | --- |
| Quadratic foundation | `quadraticNonnegative`, `lagrangeIdentity` | Prove the 2×2 diagonal/determinant criterion, including a zero diagonal, and the polynomial identity. |
| F1 | `firstPartNonnegative`, `firstPartCauchy` | The first form is nonnegative and its pairing obeys the squared Cauchy inequality. |
| F2 | `secondPartReverseCauchy` | The second form has the reverse squared inequality. |
| F3 | `complementPolarization`, `totalAtOne` | Complementary-point polarization has total \(2(p-u)^2\). |
| F4 | `weightedCauchy`, `strictWeightedCross` | Two negative vectors for the second form in the nonnegative quadrant have negative pairing. |
| Contradiction | `contradictionFromParts`, `lemmaLOrdered` | Two negative values of \(Q\) would force incompatible strict and weak square inequalities. |
| Symmetry | `cubeFormSwap`, `lemmaL`, `lemmaLMax` | Swap both parameter pairs and coordinate pairs if \(pq>ru\), giving the full statement. |

The F4 proof compares squares using the exact identities
\(LL'-T^2=AC(xy'-yx')^2\) and
\(S^2-(2Hxy)(2Hx'y')=H^2(xy'-yx')^2\), with
\(L=Ax^2+Cy^2\), \(T=Axx'+Cyy'\), and \(S=H(xy'+yx')\).
It uses strict multiplication only after deriving the necessary positivity. It therefore handles degenerate diagonals without division or unstated nonvanishing assumptions.

## Reproduction

Use the version pinned in [`formal/lean-toolchain`](../../../formal/lean-toolchain), Lean **4.24.0**:

```sh
python3 research/marcus-de-oliveira/orthogonal-cube-lean/check.py /path/to/lean-4.24.0/bin/lean
```

The checker compiles the source with `-DwarningAsError=true` into a temporary directory, requires the exact 18 public axiom reports, and rejects every transitive axiom outside `propext`, `Classical.choice`, and `Quot.sound`. It also checks the exact generic theorem type independently and compiles its rational specialization. No tracked build artifacts are written. `audit.json` records the checked source hash and declaration inventory; the checker rejects drift from that hash.

There are no `axiom`, `sorry`, `admit`, or `native_decide` proofs in the source. The arithmetic tactic `grind` produces kernel-checked proof terms. The local recursion limit of 4096 only permits elaborating those terms; it does not alter the kernel or the theorem hypotheses.
