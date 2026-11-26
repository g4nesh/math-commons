# Moving the second eigenvalue from zero to one while the fourth is unrestricted

In the fixed [four-cycle fixture](../full-four-cycle/README.md), put A=diag(x,η,y,z) and keep B fixed. For **every real η in [0,1] and all complex x,y,z**, there are explicit weights such that

    det(A+B) = Σ_p w_p(η,z) ∏_j (a[p(j)] + b[j]).

The weights work for every pair x,y simultaneously. The earlier [two entire-plane certificates](../full-four-cycle/README.md) established the endpoints η=0 and η=1. This certificate covers every point between them. The [complex-square construction](../two-middle-positive-laws/README.md) also allows a small imaginary part near zero.

## Constructing the weights

Write z=s+it. For k=0,…,3 define B_k(η)=binom(3,k)ηᵏ(1−η)³⁻ᵏ. Each entry in [certificate.json](certificate.json) gives a control index k, a zero-based permutation index p in lexicographic column-to-row order, three small integers (a,b,c), and a positive rational coefficient q. Add the following quantity to the numerator N_p for that permutation:

    B_k(η) q (a + bs + ct)².

Let H be the sum of the 24 numerators and use w_p=N_p/H. For each of the four controls, the coefficients whose square is the constant 1 sum **exactly** to one. Since the controls are nonnegative and sum to one, H≥1 everywhere on the interval, even at its endpoints or arbitrarily large z. Other squares may vanish. This establishes nonnegative weights summing to one. Our searches through degree three with the prescribed denominator 1+s²+t² did not produce a certificate. Using H allows the controls and the fourth eigenvalue to determine the normalization.

## Proving the determinant identity for all parameters

After multiplying by H, the difference between the two sides is a polynomial of degree at most four in η, total degree at most three in (s,t), and degree at most one in each of x,y. The determinant and every vertex are affine in each diagonal entry. Each numerator has degree at most three in η and two in (s,t).

The checker uses the original 4 × 4 determinant, evaluated over Gaussian rationals independently of the coefficient solver. It tests η=0,1/4,1/2,3/4,1, both x and y in {0,1}, and all ten integer points s,t≥0 with s+t≤3. Its 200 complex equalities force the polynomial to vanish identically. For the triangular grid, the four points at t=0 first force divisibility by t. At t=1 the quotient has degree at most two and three roots in s, forcing a factor t−1. Repeating at t=2 and t=3 gives zero. The five η values and the two values in each remaining variable give the other interpolation steps. This argument includes repeated spectra and collisions.

For discovery we solved linear constraints on nonnegative Bernstein control coefficients. We reconstructed 252 strictly positive coefficients over the rationals on the selected support and checked 400 real coefficient equations and four denominator-floor equations. Replaying the finished proof requires only Python's standard library.

## Reproduce or obtain weights

~~~sh
python3 -B research/marcus-de-oliveira/second-bernstein-segment/check.py
python3 -B research/marcus-de-oliveira/second-bernstein-segment/check.py 1/2 -7/3 11/5
~~~

The second command prints 24 exact fractions at η=1/2 and z=−7/3+11i/5. The eigenbasis and the opposing spectrum remain the archived ones. Extending this certificate to the whole complex η-plane remains open.
