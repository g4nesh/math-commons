# U12=0 residual inequality: a one-parameter cone certificate and its obstruction

**Status: partial.** The remaining U12=0 inequality is not proved here. This
report reduces it, at each parameter point, to the existence of one real number
m. It proves that reduction is sound by an exact identity. It shows with exact
rational witnesses that three explicit choices of m fail. No claim is made
about ★ on U12=0, on all of U(4), or about the Marcus–de Oliveira conjecture.

Input, read only: the U12=0 reduction in
[star-u12-boundary](https://github.com/g4nesh/math-commons/pull/132) (branch
`raghav/mdo-star-u12-boundary-v2`), which defines F0, Fc, F1, M and the frame.

## 1. Target

For (u, w, z) on the unit sphere with q² = w² + z² > 0, z > 0, and t > 0,
prove for every σ ≥ 0 and every unit (B, G) ∈ C²

    Φ := σ² F1 + σ Fc + F0 ≥ 0.

Write X = Re(y y*) for y = (B, G). X is a real 2×2 PSD matrix with trace 1. Then

    Fc = X : M,   F0 = t² (lᵀ X l)(X : Nn),   l = (uz, −w),
    Nn = u² q² I + t m̂ m̂ᵀ,   m̂ = (z, −uw).

Here X : A = tr(XA). `verify.py` replays the F0 factorisation. The cone of
2×2 PSD matrices is a Lorentz cone, so the S-lemma for quadratic forms on it
applies.

## 2. Cone certificate (proved sound)

Take the frame E = [l, l⊥] with l⊥ = (w, uz). Then EᵀE = l2·I with
l2 = u²z² + w², and X = E Y Eᵀ with Y ⪰ 0. Put Mp = EᵀME and Np = EᵀNnE.
Write τ = σ, Y = [[y1, y3], [y3, y2]]. The exact identity is

    Φ = [τ, y1, y3] P3 [τ, y1, y3]ᵀ + t² l2² Np22 · det Y + τ (Y : N),

    N  = [[n11, m], [m, Mp22]],
    P3 = [[F1, h1, h3], [h1, k11, k12], [h3, k12, k22]],
    h1 = (Mp11 − n11)/2,  h3 = Mp12 − m,  k = t² l2² Np.

`verify.py` checks this as a polynomial identity in all symbols.

**Sufficient condition.** Suppose Mp22 = Fc(l⊥) > 0 at a point (u, w, z, t), and
some real m and n11 ≥ m²/Mp22 make P3 ⪰ 0. Then Φ ≥ 0 at that point for all
σ ≥ 0 and all unit (B, G). Each of the three terms is nonnegative: N ⪰ 0, Y ⪰ 0,
τ ≥ 0 and Np22 ≥ 0.

The free data really is one number. Φ vanishes at τ = 0, Y = e2e2ᵀ. That forces
the det Y multiplier to be t² l2² Np22 and forces N22 = Mp22, for every
certificate of the form PSD + μ det Y + τ (Y : N). The optimal n11 for a given m
lies on n11 = m²/Mp22 whenever M is not PSD. I argued this by hand; it is not
replayed. Mp22 ≥ 0 is necessary for the target itself, since F0 = 0 at y ∝ l⊥.
It is not proved here.

## 3. Evidence that the family suffices (floats only, not a proof)

For each of 3000 random points (u, w, z, t), I chose m by a 1-D minimisation.
The points include a third with u, z ∈ [10⁻⁴, 10⁻¹] and t up to e⁸. Every point
gave F1 − hᵀK⁻¹h > 0. The minimum relative slack was 1.4·10⁻³, near u, z → 0
and t → ∞, which matches the tight limit reported in PR 132. Mp22 > 0 held on
every sample.

## 4. Obstruction: three explicit choices of m fail

Points are (u, w, z) = (2a, 1 − a² − b², 2b)/(1 + a² + b²). Both have det M < 0
and Mp22 > 0.

| point | a | b | t | passing m (normalised frame) |
|---|---|---|---|---|
| W1 | 1/491 | 1/887 | 813 | 5877/1000 |
| W2 | 1/200000 | 1/344 | 37951/10 | 17447/100 |

At both points the following choices make P3 fail to be PSD; an exact principal
minor is negative:

1. N = Fc(e2) e2e2ᵀ (m = 0, n11 = 0);
2. N rank one with h3 = 0 (m = Mp12);
3. N rank one, with m from one Newton step from m = 0 on hᵀ adj(K) h.

At each point the listed rational m (multiplied by l2 in the unnormalised frame),
with n11 = m²/Mp22, makes P3 positive definite exactly. So the certificate family
is non-empty there, but m must depend on the point in a nonlinear way. The float
optimum is sensitive: at W2, m = 174 fails and m = 174.47 passes.

## 5. Failed approaches

* N = Fc(e2)e2e2ᵀ: fails on about 90% of float samples, including some with
  det M < 0 (W2).
* h3 = 0 with minimal n11 (the "complete M along l⊥" choice): fails on most
  samples with det M < 0. Rank-one κ completion from PR 132 is a sub-case of this
  family with μ = 0.
* One Newton step for m: overshoots at generic points, for example near
  (0.23, 0.79, 0.57), t = 40.

## 6. Next step

The optimal m solves the cubic ∂m (hᵀ adj(K) h) = 0 with h1 = (Mp11 − m²/Mp22)/2.
Two routes remain. (a) Eliminate m: "∃ m with P3 ⪰ 0" is a semialgebraic
condition in (u, w, z, t). A resultant or discriminant form of it could be
certified by Bernstein coefficients in (U, W, Z, t) on a subdivision. (b) Use a
piecewise rational selector, one per region: |u|, |z| small with t large, and
the rest. The adapted coordinates t → 1/t near the tight corner suggested in
the task fit here. Mp22 > 0 also still needs an exact proof.

## 7. Reproduce

From the repository root, with Python 3.11 and SymPy 1.14:

```sh
python -I -B research/marcus-de-oliveira/star-u12-residual/verify.py
```

It prints 20 `ok` lines and `ALL CHECKS PASSED`, then exits 0. It takes about
30 s. Every decision is an exact polynomial identity or an exact rational sign.
The float searches in §§3–4 only chose candidates and are not in the script.

Derived by claude-opus-5-5-residual (Claude Opus 5.5, operator Raghavsk24).
Same-owner checks are not independent review.
