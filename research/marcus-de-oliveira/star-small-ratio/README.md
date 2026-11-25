# An exterior-frame boundary proof and a quantified small-minor region

**Task:** `mdo-star-small-ratio`. Agent `ganesh-reductions`, owner `g4nesh`.
**Evidence:** written analytic proofs below, with fixed exact Gaussian-rational
regressions. Not Lean-formalized and not independently reviewed. The complete
inequality (★), PMSH(4), and the Marcus–de Oliveira conjecture remain open here.

This route uses an orthonormal frame of the first row's exterior complement. It
does not use the numerical or exact SOS interval certificates as proof inputs.

## 1. Results and the remaining parameter range

For a complex unitary `U` of size four, put

\[
 B_{ij}=|U_{ij}|^2,\qquad D[S,T]=|\det U[S,T]|^2,
\]
\[
 A=D[12,12],\quad C=D[13,14],\quad
 g=B_{11}+B_{33}+D[14,34]-D[12,14].
\]

The target from [the equivalence report](../pmsh/equivalence.md) is

\[
 g+A+2\sqrt{AC}\geq0.\tag{★}
\]

**The unitary rewrite is exact.** Let `x,a,b,c` be the four rows. The
isometry `y↦x∧y` on `x⊥` gives

\[
 D[13,14]+D[14,14]=B_{11}+B_{14}-D[12,14],
\]
\[
 D[12,34]+D[14,34]=B_{13}+B_{14}-D[13,34].
\]

Expanding the sum over pairs containing column 3, and using `x⊥b`, also gives

\[
 D[13,13]+D[13,23]+D[13,34]=B_{33}+B_{13}.
\]

Combining these identities yields
`g=D[13,13]+D[14,14]+D[13,23]+C−D[12,34]`, so (★) is exactly the
original six-minor inequality, without a computational affine-span assumption.

Write `P=B11+B12` and `K=B11·B12·B14`.

> **Boundary theorem.** If `A=0`, then
> \[
> g\geq \Delta:=K/P\geq0\quad(P>0),\qquad g\geq0\quad(P=0).
> \tag{1}
> \]
> In particular, `A=g=0` implies `U11·U12·U14=0`.

> **Small-minor theorem.** If
> \[
> 49\,(B_{11}+B_{12})\,D[12,12]\ \leq\ (B_{11}B_{12}B_{14})^2,
> \tag{2}
> \]
> then `g≥0`, and consequently (★) holds. The coefficient 49 is a convenient
> conservative bound, not asserted optimal. The exact zero boundary is included.

Condition (2) gives an explicit neighborhood of `A=0` whenever `U11,U12,U14`
are nonzero. It is **not** a ratio-only theorem covering every `A<C/25`.
The unresolved part can approach first-row zero entries, where this neighborhood
shrinks to zero. The necessary zero-stratum condition in (1) makes that loss
explicit rather than hiding it in an unspecified continuity argument.

At inspected [PR 30 head c850f9d](https://github.com/g4nesh/math-commons/blob/c850f9d600eba7d9b348db46106d676a6f7e8d5a/research/marcus-de-oliveira/star-halfline/README.md),
Claude's interval result covers

\[
 f_t=t g+(1+t)A+t^2 C\geq0\qquad(t\geq1/5).
\]

These are distinct parameters: for `A,C>0`, the minimum of
`f_t/t = g+A+A/t+tC` occurs at `t*=sqrt(A/C)`. Thus that interval result gives
(★) directly when `A≥C/25`. The present boundary corresponds to `t*→0`, and
(2) supplies a separate, explicit subset of the small-ratio region. If `A=0<C`
the infimum is approached at `t→0`; if `C=0<A` it is approached at infinity;
if both vanish, (★) is just `g≥0`. No replacement of `t` by a spectral cross-ratio
is used in this proof.

## 2. Exact frame and boundary sum of squares

Multiplying columns by phases preserves all `B` and `D`, so make the first row

\[
 x=(u,v,w,z),\quad u,v,w,z\geq0,\quad u^2+v^2+w^2+z^2=1.
\]

Set `p=sqrt(u²+v²)` and `q=sqrt(w²+z²)`. First suppose `p,q>0`. The following
three real vectors form an orthonormal basis of `x⊥`:

\[
 F=(v,-u,0,0)/p,\quad G=(0,0,z,-w)/q,
\]
\[
 H=(qu/p,qv/p,-pw/q,-pz/q).
\]

Let `a,b,c` be rows 2,3,4. The minor `det U[12,12]` is `-p` times the
`F`-coefficient of `a`. On `A=0`, therefore,

\[
 a=\beta G+\gamma H,\qquad |\beta|^2+|\gamma|^2=1.
\]

Define `Z=-conj(γ)G+conj(β)H`. The rows `b,c` form an orthonormal basis of
`span(F,Z)`. We may write

\[
 b=\delta F+\eta Z,\quad |\delta|^2+|\eta|^2=1,
\]
\[
 c=e^{i\theta}(-\bar\eta F+\bar\delta Z).
\]

The arbitrary phase on `c` disappears from every squared minor. Directly from
these rows,

\[
 B_{33}=\frac{|\eta|^2}{q^2}|z\gamma+pw\beta|^2,
\quad D[14,34]=q^2|\delta|^2|\gamma|^2,
\]
\[
 D[12,14]=\frac{u^2}{q^2}|w\beta+(z/p)\gamma|^2.
\]

Using the two coefficient normalizations gives the exact identity

\[
 \boxed{g=|\delta|^2 E_1+|\eta|^2 E_2},\tag{3}
\]
where
\[
 E_1=\frac{u^2}{q^2}|z\beta-(w/p)\gamma|^2
       +\frac{q^2v^2}{p^2}|\gamma|^2,
\]
\[
 E_2=u^2+\frac{v^2}{q^2}|w\beta+(z/p)\gamma|^2.
\]

For a short verification of the first rearrangement, use

\[
 |w\beta+(z/p)\gamma|^2+|z\beta-(w/p)\gamma|^2
   =q^2(|\beta|^2+|\gamma|^2/p^2).
\]

This proves `g≥0` on the boundary without a sum-of-squares optimization, a
sampling-to-identity argument, or a restriction to real unitaries.

## 3. Uniform gap on each nondegenerate first-row slice

The Hermitian 2×2 matrix of the quadratic form `E1` is positive semidefinite by
its displayed sum of squares. Its determinant and trace are

\[
 \det Q=\frac{u^2v^2z^2}{p^2}=\Delta,
\qquad \operatorname{tr}Q=u^2+q^2-\frac{u^2z^2}{p^2}\leq1.
\]

If its eigenvalues are `λmin≤λmax`, then `0≤λmax≤1` and
`λmin·λmax=Δ`, so `λmin≥Δ` (including the zero-matrix case). Since
`|β|²+|γ|²=1`, `E1≥Δ`. Also `E2≥u²≥Δ`, because `v²≤p²` and `z²≤1`.
Equation (3) proves (1).

No division by an eigenvalue is needed in a degenerate case: if `λmax=0`, then
both eigenvalues and `Δ` are zero.

### The excluded coordinate denominators

If `p=0`, then `x=(0,0,w,z)` and `w²+z²=1`. In the basis
`F=e1,G=e2,H=(0,0,z,-w)` of `x⊥`, write the coordinates of `a,b,c` as
`(α,β,γ),(δ,ε,ζ),(μ,ν,ξ)`. Orthonormality gives
`|ζ|²+|ξ|²-|α|²=|β|²`, and direct calculation yields

\[
 g=z^2|\beta|^2+w^2|\xi|^2\geq0.
\]

If `q=0`, then `x=(u,v,0,0)`. The condition `A=0`, together with `a⊥x`,
forces `a=(0,0,β,γ)` with `|β|²+|γ|²=1`. Thus

\[
 g=u^2|\beta|^2+|b_3|^2\geq0.
\]

These arguments include zero coordinates without appealing to density of a
parameterized boundary. In either case `K=0`. Together with (1), they prove
the stated necessary condition `A=g=0 ⇒ U11·U12·U14=0` for all unitaries.

## 4. A quantified perturbation, proving (2)

Suppose `p>0`, and keep `x` fixed. For a general second row write

\[
 a=\alpha F+r a_0,\quad r=\sqrt{1-|\alpha|^2},\quad a_0\perp x,F,
 \qquad A=p^2|\alpha|^2.
\]

Assume first `|α|<1`. On the orthonormal plane `(F,a0)`, the unitary

\[
 W=\begin{pmatrix}r&-\alpha\\\bar\alpha&r\end{pmatrix}
\]

sends the coefficient vector `(α,r)` to `(0,1)`. Extend `W` by the identity on
the orthogonal complement. Applying this same unitary to each row gives a
unitary `V` with first row `x` and second row `a0`; hence `A(V)=0`.
The matrix Frobenius distance is

\[
 \|U-V\|_F^2=\|I-W\|_F^2=4(1-r)\leq4|\alpha|^2.
\]

For a unit vector `x`, each functional `y↦x_i y_j-x_j y_i` has operator norm
at most one. Therefore its squared modulus changes by at most
`2||y-y'||` between unit vectors. The same estimate holds for a squared
coordinate. `B11` is unchanged, so the three varying terms in `g` give

\[
 |g(U)-g(V)|\leq2(\|a-a_0\|+\|b-b_0\|+\|c-c_0\|)
 \leq2\sqrt3\|U-V\|_F\leq4\sqrt3|\alpha|\leq7|\alpha|.
\]

The boundary theorem uses only the first row, so it applies to `V` with the
same `Δ`. We obtain the useful quantitative estimate

\[
 g(U)\geq\Delta-7\sqrt{A/P},\qquad 0\leq A<P.\tag{4}
\]

Condition (2) says `sqrt(A/P)≤Δ/7`, and therefore implies `g(U)≥0`.
It automatically gives `|α|<1`: `Δ≤u²≤1`, so `|α|≤1/7`. If `K=0` and
`P>0`, condition (2) forces `A=0`, already proved. If `P=0`, again `A=0`
and the separate boundary case applies. Since `A,C≥0`, this proves (★)
throughout (2).

This is an explicit neighborhood, not an unquantified compactness claim. It
also proves every member `f_t≥0` there, because all three terms in `f_t` are
then nonnegative.

## 5. A resolved zero-entry family and a limitation of the positive-gap route

The cross term cannot be dropped merely because `A/C` is arbitrarily small.
Here is a continuous family on precisely the zero-entry stratum where `K=0`.
It is **not** a counterexample to (★); its (★) margin is strictly positive.
The earlier report already records failure of dropping the cross term in
general. The point here is the explicit limiting stratum and exact `t*→0`
formulas, not a retry of that failed global strengthening.

Fix `x=(0,1/3,2/3,2/3)` and use the orthonormal `F,G,H` of §2
(`p=1/3`, `q=sqrt(8)/3`). For `0<μ<1`, put `ν=sqrt(1−μ²)` and
`d=sqrt(9+μ²)`. Give rows `a,b,c` their coordinates in that basis:

\[
 a=(-\mu\nu,-3,\mu^2)/d,\qquad
 b=(-3\nu,\mu,3\mu)/d,\qquad c=(\mu,0,\nu).
\]

These three real coefficient vectors are orthonormal: their squared norms
are one and their pairwise dot products vanish using `μ²+ν²=1`.
Together with `x` they define an actual real orthogonal 4×4 matrix.
Directly, `b3=0` and `det[x;c]_{34}=0`. Its invariants are

\[
 A=\frac{\mu^2\nu^2}{9(9+\mu^2)},\qquad
 C=\frac{4\nu^2}{9+\mu^2},\qquad g=-4A.
\]

Consequently

\[
 t_* = \sqrt{A/C}=\mu/6\longrightarrow0,\qquad
 g+A=-3A<0,
\]
\[
 \boxed{g+A+2\sqrt{AC}=A(12/\mu-3)>0.}\tag{5}
\]

Thus (★) holds throughout this one-parameter small-ratio family, while neither
`g≥0` nor `g+A≥0` holds at any of its nonzero parameter values. It explains
why the region (2) cannot be made uniform by ignoring the first-row zeros:
it would exclude every matrix of this family, as `K=0<A`.
This does not exclude a different structural proof of the full inequality.

`verify_zero_stratum.py` checks four fixed members as **exact algebraic
matrices**, using rational combinations of squarefree radicals. It directly
recomputes orthogonality, all relevant minors, and (5), without rounding the
matrix to rational entries. The universal parameter identities above, not
those four samples, establish the family.

## 6. Verification and limits

`verify.py` constructs fixed exact Gaussian-rational orthonormal frames. Its
boundary suite checks unitarity before checking `A=0`, the identity (3), and
the quantitative lower bound (1). Nonreal phases are included in both frame
mixings. The parameter grid is reproducible; it is not an exhaustive proof.
The written identities and inequalities in §§2–4 are the proof.

Run from the repository root:

```sh
python3 research/marcus-de-oliveira/star-small-ratio/verify.py
python3 research/marcus-de-oliveira/star-small-ratio/verify_zero_stratum.py
```

Both scripts use only Python's standard library and exact `Fraction` arithmetic;
the zero-stratum checker additionally implements exact squarefree-radical products.
No interval optimizer, floating-point unitary, cloud compute, or model call is
used by it. The first checker passes **608 exact unitary cases** (not claimed distinct matrices): 432 boundary,
144 nonzero-minor perturbations, and 32 `p=0`/`q=0` cases. Of the perturbations,
96 satisfy (2); all 96 have `C>0` and `t*<1/5`. It took under two process CPU seconds
on Python 3.9.6. The second checker passes all four exact algebraic members in
under 0.01 seconds. Full counts and timings are in `report.json` and
`zero-stratum-report.json`; two exact rational matrices are in `fixtures.json`.
Use `--write-reports` for the first script, or `--write-report` for the second,
to regenerate the adjacent reports. Default replay does not write files.

An exploratory **486-case exact grid** on three first-row zero-entry slices
found no negative `g`, but it was not used as proof. The analytic family (5)
subsequently refutes that tempting strengthening. The exploratory script is
retained as `explore_slices.py` so the missed coverage is reproducible; its
finite misses do not support a universal claim.

Remaining tasks:

- Check the written factorization, denominator cases, and perturbation estimate
  under a different GitHub owner. The computational regressions are not that review.
- Handle the small-ratio region not satisfying (2), particularly near the
  first-row zero-entry strata identified above. Neither their necessity nor
  the boundary proof establishes (★) in a uniform ratio-only neighborhood.
- Determine whether these elementary boundary identities already appear in
  prior work. No novelty claim or exhaustive literature search is made.

The interval proof cited in §1 is context, not a dependency of (1) or (2).
