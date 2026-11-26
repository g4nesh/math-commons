# Two complex eigenvalues and one whole real line

For the archived real orthogonal four-cycle

$$
U=\frac1{625}\begin{pmatrix}
-31&492&-240&300\\
-300&225&500&0\\
240&-180&225&500\\
-492&-256&-180&225
\end{pmatrix},\qquad b=(3-2i,2,1+2i,-i),
$$

put $B=U\operatorname{diag}(b)U^*$ and $A=\operatorname{diag}(x,t,y,0)$.
**For every real $t$, one explicit probability distribution on permutations works simultaneously for every $(x,y)\in\mathbb C^2$.** Consequently the original MDO inclusion holds for this entire five-real-dimensional family. The weights depend on $t$, but not on $x$ or $y$.

The requested interval atlas simplifies to one global formula: its numerator weights are sums of positive multiples of real squares, and its denominator never vanishes on the real projective line. This includes $t=-2$ and $t=0$, where some coefficient bases used during discovery lose rank. No limiting argument at those two finite points is needed.

## The global formula

A permutation $p$ maps column to row; define $V_p(x,t,y)=\prod_{j=0}^3(a_{p(j)}+b_j)$, with $a=(x,t,y,0)$. The machine-readable [certificate](certificate.json) gives 29 triples $(p,r,c)$, all $c\in\mathbb Q_{>0}$ and $r\in\mathbb Q$. Put

$$q_p(t)=\sum_{(p,r,c)}c(t-r)^2,\qquad w_p(t)=\frac{q_p(t)}{Q(t)},$$

where an absent permutation has numerator zero. Exactly 14 of the 24 permutations are used. The sum of their numerators is

$$
Q(t)=t^2+\frac{67729206836}{79985617449}t
          +\frac{721641819716}{79985617449}.
$$

Every $q_p$ is nonnegative on the entire real line. For $Q(t)=a+bt+ct^2$,

$$c=1,\qquad 4ac-b^2=
\frac{226296620649383138167040}{6397698998697773267601}>0.$$

Completing the square proves $Q(t)>0$ for every real $t$. Thus the $w_p$ are nonnegative and sum to one. The decisive exact polynomial identity is

$$
\sum_p q_p(t)V_p(x,t,y)
=Q(t)\det\bigl(\operatorname{diag}(x,t,y,0)+B\bigr).
$$

The verifier constructs both sides as polynomials over Gaussian rationals in the three formal variables $x,t,y$ and checks every coefficient. Their coordinatewise degrees are at most $(1,3,1)$: 16 complex coefficient slots, or 32 real slots. This is an identity over $\mathbb C^3$, while **nonnegativity of these weights is asserted only for real $t$**. Dividing by positive $Q(t)$ proves the stated common-weight theorem for arbitrary complex $x,y$.

## The chart at infinity

For a real projective point $[u:v]$, define

$$q_p^h(u,v)=\sum_{(p,r,c)}c(u-rv)^2,\qquad
Q^h(u,v)=u^2+\frac{67729206836}{79985617449}uv
                  +\frac{721641819716}{79985617449}v^2.$$

The same discriminant proves that $Q^h$ is positive for every nonzero real pair $(u,v)$. The weights $q_p^h/Q^h$ therefore extend continuously over all of $\mathbb{RP}^1$ and are unchanged by rescaling $(u,v)$. On $v\ne0$ they give the finite formula at $t=u/v$. The reciprocal chart uses $(u,v)=(1,z)$, so $z=0$ is covered directly.

At $[1:0]$, the weights are the sums of the positive constants $c$ assigned to each permutation and sum to one. Homogenizing the displayed polynomial identity, with determinant and vertex pencils of degree one in $t$, proves the corresponding equality for their leading coefficients. Infinity is a coefficient-pencil limit, **not a finite matrix with an infinite eigenvalue**.

## Proof replay and discovery

Run from the repository root with Python 3 and no third-party package:

```sh
python3 -B research/marcus-de-oliveira/two-free-real-atlas/check.py
python3 -B -m unittest discover -s research/marcus-de-oliveira/two-free-real-atlas -p 'test_*.py' -v
python3 -B research/marcus-de-oliveira/two-free-real-atlas/check.py --weights -2
python3 -B research/marcus-de-oliveira/two-free-real-atlas/check.py --infinity
```

The standalone checker verifies the displayed $U$ is exactly orthogonal, constructs $B$ directly, expands the original $4\times4$ determinant by its 24 signed terms, and independently constructs all pairing polynomials. It checks positive square masses, the numerator sum, the exact positive discriminant and the full polynomial identity. No optimizer, polynomial library, external source file, sampled positivity or floating-point tolerance is used in the proof replay. The seven regression tests additionally compare with the archived Gaussian-elimination determinant implementation at non-Boolean complex spectra, exercise both singular discovery points and infinity, and reject corrupted certificates and inexact API inputs.

`weights_at(t)` accepts an integer or `Fraction`. `weights_homogeneous(u,v)` accepts a nonzero pair of such exact reals and includes infinity at `(1,0)`. These API inputs are exact rational representatives; the polynomial and positivity proof covers all real parameters, including irrational ones.

Discovery first used linear programming to find feasible bases at selected parameters; it detected coefficient rank drops at $-2$ and $0$. Instead of continuing these bases across intervals, a second linear program searched positive combinations of the primitives $(t-r)^2$, with candidate centers $r\in\{-10,-3,-2,-1,-1/2,0,1/2,1,2,3,10\}$ and an optional constant primitive. The resulting support used 29 squares. Solving its coefficient equations over $\mathbb Q$ gave the stored positive masses; the standalone expansion above independently rechecks their identities. Numerical discovery supplies no premise of the proof.

## Relationship to the uniform obstruction

The archived [fourth-eigenvalue boundary](../fourth-boundary/README.md) excludes one fixed law valid for all three complex eigenvalues at $a_4=0$, using an exact strict separator. The present law changes with $t$. It therefore establishes the original determinant inclusion on a large family despite that stronger uniform obstruction; it does not contradict the separator. Even a law constant for all real $t$ and all complex $x,y$ would imply the prohibited three-variable coefficient identity, so dependence on $t$ is essential.

This is a fixed-$U$, fixed-$b$, fixed-$a_4$ theorem. It does not cover arbitrary complex $t$, all normal pairs, or the general MDO conjecture. There is no claim of literature novelty, Lean formalization, or independent-owner review. The next distinct question is whether direction-dependent formulas for $t=(1+ik)u$ and $t=(k+i)u$, with real $u,k$, can cover every complex $t$ by a finite exact atlas in the direction parameter.
