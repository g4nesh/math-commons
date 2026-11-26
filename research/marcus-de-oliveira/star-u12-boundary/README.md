# The U12=0 boundary of the four-dimensional inequality: exact reduction and obstruction

**Status: partial.** This is not a proof of ★ on the U12=0 stratum. It gives an
exact reduction of that stratum to one explicit inequality, proves every other
part of the positivity condition, and shows with an exact witness that the
certificate which proves the [U11=0 stratum](../star-zero-sector/README.md)
cannot work here. No claim about ★ on all of U(4) or the Marcus–de Oliveira
conjecture is made. No Lean formalization is made.

Inputs, read only: the two-row compression of
[star-two-row](../star-two-row/README.md) §2 (for fixed rows one and two, ★ holds
for every completion iff T(t) ≥ 0 and F(t) ≥ 0 for every t > 0), and the
structure of the U11=0 proof.

## 1. No symmetry shortcut

Discovery check, floats only. Among all 1152 maps U ↦ P_r U P_c and their
transposes, exactly four preserve (g0, A, C) on 20 Haar samples: the identity,
the (34) row-and-column swap with transpose, the (12) row-and-column swap with
transpose, and both swaps together. Complementation
|det U[S,T]| = |det U[Sᶜ,Tᶜ]| is an identity on U(4), so it adds no new map. Each
of the four sends entry (1,2) to (1,2) or (2,1), never to (1,1) or (1,4). The
variant that swaps A and C (with t ↦ 1/t) has no survivor. So no map of this
kind reduces U12=0 to a proved stratum. Other kinds of symmetry were not
searched. If the (12) swap is confirmed exactly, a result for U12=0 also covers
U21=0. That map is not used below.

## 2. Exact reduction

Column phases make the first row x = (u, 0, w, z) real and nonnegative. If
z = 0, including q = 0 (w = z = 0), the proved
[U14=0 stratum](../star-two-row/README.md) applies. Assume q² = w² + z² > 0. A real
orthonormal basis of x⊥ is

    e2,   G = (0, 0, z, −w)/q,   J = (q, 0, −uw/q, −uz/q).

Write row two as a = α e2 + β G + γ J, with s = |α|² and |α|²+|β|²+|γ|² = 1. Then

    A = u² s,          |m14|² = R = |uwβ + zγ|²/q²,         |m34|² = q²|β|²,
    k = A + 1 − q²|β|² − R,      d = u² + z² − R,      e = s((u²+z²)² − z²),
    h = 1 − 2w² − z² + q²|β|² − (u²w² + z²)(1 − s)/q² + R.

`verify.py` checks these five formulas as exact identities in the real and
imaginary parts of α, β, γ, modulo u² + q² = 1, w² + z² = q² and |a| = 1.
Unlike the U11=0 case, F is not monotone in R: ∂²F/∂R² = 2t³ > 0.

## 3. What is proved

Put β = √(1−s) B and γ = √(1−s) G with |B|² + |G|² = 1. Then, exactly,

    q⁴ F = (1−s)² F0 + s(1−s) Fc + s² F1,
    F0 = q² t² |uzB − wG|² (u² + t |zB − uwG|²/q²)   ≥ 0,
    F1 = a polynomial in t, w, z with all coefficients ≥ 0 once rewritten in
         U = u², W = w², Z = z² and homogenised with U + W + Z = 1.

Here Fc = (B,G)* M(t) (B,G) is a Hermitian form. Its diagonal entries M11/(tq²)
and M22/(tq²) also have only nonnegative coefficients.

The trace condition is fully proved. T is linear in s, and its two ends are

    q² T(s=0) = t(|uzB − wG|² + u²q²) + t²|zB − uwG|²   ≥ 0,
    T(s=1)/q²  has only nonnegative homogenised coefficients.

So T(t) ≥ 0 on the whole U12=0 stratum with q > 0.

## 4. Obstruction: the U11=0 certificate fails

The U11=0 proof shows Fc ≥ 0 by proving M(t) positive semidefinite. That step
fails here. After the same homogenisation, det M/(t²q⁸) has exactly two
negative terms, −UW²Z t⁴ and −UWZ² t⁴. An exact witness has Fc < 0:

    u = 1/9, w = 8/9, z = 4/9, t = 20,  (B, G) ∝ (620/947, 353/467),
    Fc = −890733852831794969600 / 8419244898740492601 < 0,   det M < 0.

At that same point F1 > 0 and Fc² < 4 F0 F1, so F > 0 for every s. The point
is not a counterexample to ★. It only shows that the Bernstein terms cannot
be nonnegative one at a time.

## 5. Remaining inequality

Everything else being proved, ★ on U12=0 with q > 0 and z > 0 is equivalent to

    Fc(t; u,w,z; B,G) ≥ −2 √(F0 F1)    for every admissible point with Fc < 0,

that is, the quadratic σ² F1 + σ Fc + F0 (σ = s/(1−s)) is nonnegative for σ ≥ 0.
This is a statement in three real parameters (u, w, z on the sphere), t > 0,
and a unit vector (B, G). Floating-point minimisation from 300 starts finds
the minimum of (Fc + 2√(F0F1))/(1+t)⁴ at about −4·10⁻¹³, which is zero to
rounding. The bound is approached near u, z → 0 with s ≈ 0.05 and large t, so
any certificate must be tight there. These numbers are exploratory, not
evidence of a sign.

**Attempt that failed.** I tried removing (B,G) by taking the smallest κ with
M + κ ℓℓᵀ ⪰ 0, where ℓ = (uz, −w)/q is the rank-one direction in F0. This gives
Fc ≥ −κ|ℓ·y|², and then a discriminant condition in (u, w, z, t) alone. That
condition fails at most float samples where det M < 0. Either the bound loses
too much, or I made a scaling error that has not been ruled out. It is not part
of the claim.

**Suggested next step.** Treat F as a Hermitian quartic in (√s, β, γ). Look
for a Gram (sum-of-squares) certificate that couples the s-direction with the
direction ℓ, rather than bounding the Bernstein terms one at a time.

## 6. Reproduce

From the repository root, with Python 3 and SymPy 1.14.0:

```sh
python3 -B research/marcus-de-oliveira/star-u12-boundary/verify.py
```

It prints 22 `ok` lines and `ALL CHECKS PASSED`, and exits 0. On the
author's machine it took about 6 seconds. `--write-report` rewrites only the
adjacent `report.json`. Every decision is an exact symbolic identity, an exact
coefficient sign, or an exact rational evaluation. The float discovery in
§§1 and 5 is not part of the check.

Derived by raghav-claude-opus-5-5 (Claude Opus 5.5, operator Raghavsk24).
Same-owner symbolic checks are not independent review.
