# T4: a six-gap scalar certificate and explicit margin

Task `mdo-zheng-t4-chart-sign`, agent `ganesh-reductions`, owner `g4nesh`.
Written proof and an independent exact symbolic replay; submitted for another
owner's review. This is a scoped audit/refinement of Zheng's T4 argument, not a
new claim of priority or a proof of all contact cases.

The main result below needs only **six** nonnegative gaps in the distinct,
nonreal chart, and gives a quantitative positive coefficient margin. It verifies
the scalar input consumed by the separate SO(4) four-ray certificate.

## Current sources and what changed during this audit

The original source is Zheng's five-page
[T4 note](https://github.com/DongzheZheng/The-Marcus-de-Oliveira-conjecture-for-n-4/blob/53ec16567d2307266ce52a8e9678410b2457cc27/proofs/manin_transfer_20261021/root_moc4_T4_uniform_20261022.pdf),
pp. 3–4. Its formulas and proof idea are credited to that note. No author
program or certificate is executed or imported by this replay.

During the interruption of this task, PR66 supplied
[first](../zheng-n4-audit/docs/root_moc4_T4_uniform_20261022.md) and
[second](../zheng-n4-audit/docs2/root_moc4_T4_uniform_20261022.md) T4 reports.
They independently obtain the denominator factors and the necessary `Q != 0`
identity below, and correctly separate positivity of an arrowhead matrix from
its unproved-in-the-note identification with an SO(4) response. These are now
corroborated results, not outstanding gaps or discoveries of this task.
The [SO(4) atoms](../zheng-n4-audit/foundations/so4-atoms.md) and
[three-plus-three-circle proof](../zheng-n4-audit/foundations/three-plus-three-circle.md)
provide additional dependencies. Claude's merged
[unitary factorization](../unitary-factorization/README.md) supplies a different
route from real orthogonal to unitary squared-minor statements. None of those
larger proofs is silently certified by the scalar replay here.

The refinements here are the explicit six-gap normalized hypothesis, an exact slack
identity, and lower bounds that make the strictness inspectable. The script
starts from the 24 matching products rather than from a saved certificate.

## Domain and conclusion

Indices of permutations of `(0,1,2,3)` are lexicographic, beginning at zero.
Normalize the eight distinct projective points to

`A=(z,0,1,x), C=(w,v,infinity,t*v)`.

For a permutation `p`, set

`F_p = - product_{i: p(i)!=2}(A_i-C_{p(i)})`,
`R_p=(F_p-F_0)/(v*(z-w))`.

Assume the six support contacts `T4={0,1,4,5,7,8}` for a nonzero support
multiplier. The elementary contact argument below makes `x,t,r,u` real. We
consider `Im(v)!=0`. Rescale so `rho=mu*Im(v)=-1`, and define real gaps

`g_p=(R_p(v)-R_p(conj(v)))/(v-conj(v))`.

**Six-gap theorem.** In this chart, assume only

`g_j >= 0` for `j in {6,10,14,15,16,19}`, and `d=g_2-g_13>0`.

Put `a=g_15, b=g_19, c=g_10`. Then

`a*b*(c-d) - d^2*(a+b+c+d) > 0`.

The normalization is legitimate for an original supported source in this
branch: before rescaling the displayed divided differences are multiplied by
`-rho`; `a+d=(x-1)c>0` gives `c>0,x>1`, and
`g_2=-rho*x/|v|^2>=0` forces `rho<0`. The theorem itself is stated on the
normalized chart, so it does not assume full support or claim an equivalence
between six gaps and all 24 gaps. The original source needs the additional
condition `g_2>=0` only to establish the normalization sign. Once normalized,
`x>1` already makes `g_2=x/|v|^2>0`.

## Chart, signs and nonvanishing

The contact ratios are `0,1-x,t-x,x*(t-1),r-x+1,u`. In particular,
`R_5-R_4=-t*(R_1-R_0)`. Both corresponding support differences are purely
imaginary and nonzero, so `t` is real. Since `R_4-R_1=t-1!=0`, the common
multiplier `kappa*v*(z-w)=i*mu` is purely imaginary with real nonzero `mu`.
Then `x,r,u` are real. Distinctness supplies each nonzero factor used here.

The original definitions and their inverse are

```
r = z*(w-v)*(1-t*v)/(v*(z-w))
u = -(x-t*v)*(z*(1+v-w)-v)/(v*(z-w))
Dz = r*t*v^2-r*v*x+t^2*v^2-t*u*v^2-t*v*x-t*v+u*v+x
Dw = -r*t*v+r*x+t*u*v-u
z = r*v*(t*v-x)/Dz
w = (x-t*v)*(r-t*v+1)/Dw.
```

Substitution in the original variables gives

```
Dz = -(v-w)*(x-t*v)*(t*v-1)/(z-w)
Dw = -(z-v)*(x-t*v)*(t*v-1)/(v*(z-w)).
```

Thus neither denominator vanishes for distinct points. This rules out a lost
rank-drop branch; merely observing that a formula has a denominator would not.

Write `y=x-1, h=1-t, H=h+u, L=y-r+u`, and

```
E = H*(1-t*x)-h*r*y
J = H*(1-t*x)+h*y
k = g_16, e = g_14
K = t*y*b-x*h*d.
```

The raw-product identities are

```
g_2 = x/|v|^2
g_6 = r*t*y/|t*v-1|^2
c = -t*u*x/|t*v-x|^2
a = y*c-d
k = h*c-t*(b+d)
E*e = t*y*H*b-J*d
J-E = h*y*(r+1),   J-x*h*H = -u*y.
```

As `a>=0,d>0,c>=0`, we have `c>0,y>0`. Hence `r*t>0` and `t*u<0`;
`r=0` would give a collision from its original definition. If `t<0`, then
`r<0,u>0,h>0,H>0,E>0,J>0`, so `E*e=t*y*H*b-J*d` contradicts `e,b>=0`.
Consequently `t>0,r>0,u<0`. The identity for `k>=0` forces `h>0`.
Thus `0<t<1`, with no extra sign hypotheses imported from a numerical search.

## Proof that K is positive

If `E>=0`, then `J=E+h*y*(r+1)>0`. The second source identity forces `H,b>0`
and `t*y*H*b>=J*d`. Therefore

`K >= (-u)*y*d/H > 0`.

This includes `E=0`; no division by `E` occurs.

For `E<0`, define

`Q(v)=-t*(r+t-u)*v^2+(r*t+r*x+t*x+t-u)*v-x*(r+1)`.

Crucially `Q(v)=Dz*(z-v)/v!=0`. The exact replay reconstructs real symmetric
2-by-2 matrices `MD,MK` with

```
d*|Q(v)|^2 = (1,conj(v))*MD*(1,v)^T
K*|Q(v)|^2 = (1,conj(v))*MK*(1,v)^T
det(MD) = -r^3*t^2*u*x*y*L*E
det(MK) = -r^2*t^2*u^2*x^2*h*y*L*E
(1,x/t)*MD*(1,x/t)^T = r*u*x*y*(h*r-H*y)
(1,x/t)*MK*(1,x/t)^T = r*u^2*x^2*y*(y-h).
```

Suppose `L<=0`. If `h*r-H*y<=0`, then `H>0` and `r<=H*y/h`.
The identity `E-H*(h-y)*(1+y)=h*y*(H*y/h-r)>=0` forces `y>h`.
But `L=(-u)*(y-h)/h+(H*y/h-r)>0`, a contradiction. Hence `h*r-H*y>0`.
Now `det(MD)>=0`, while its value on `(1,x/t)` is strictly negative.
A real symmetric 2-by-2 matrix with these properties is negative semidefinite,
including determinant zero. This contradicts `d*|Q(v)|^2>0`. Thus `L>0`.

If `y<=h`, then `L>0` gives `H=h+u>h+r-y>0`, and
`E=(h-y)*(H+h*y)+h*y*L>0`, again a contradiction. Therefore `y>h`.
Now `det(MK)>0` and its value on `(1,x/t)` is strictly positive, so `MK` is
positive definite. This proves `K>0` for complex `v` too.

It also gives the explicit bound

`K >= det(MK)*(1+|v|^2)/(trace(MK)*|Q(v)|^2) > 0`.

Indeed its smaller eigenvalue is `det(MK)/lambda_max >= det(MK)/trace(MK)`.
This is an inequality over the entire parameter region, not sample evidence.

## Exact slack and boundaries

Using `a=y*c-d` and `k=h*c-t*(b+d)`, the capacity slack equals exactly

`a*b*(c-d)-d^2*(a+b+c+d) = (c/h)*((b+d)*K+y*b*k)`.

Every factor needed on the right is now positive or nonnegative, and
`c*(b+d)*K/h>0`. Combining with either explicit bound for `K` gives a strict
margin without treating the middle coefficients as independently positive.

The proof covers `E=0`, `k=0`, `e=0`, and singular `MD` in the contradiction.
It does not extend this strict assertion to `d=0`, `Im(v)=0`, or collisions.
Those are separate branches of the full T4 proof. Positivity of the displayed
arrowhead pencil alone is **not** an identification with the determinant
response; the SO(4) atom certificate remains a separate input.

## Reproduce

Python 3.9.6, SymPy 1.14.0 (python-flint arithmetic backend in this run).
From the repository root:

```sh
python3 research/marcus-de-oliveira/zheng-t4-chart-sign/check.py
python3 -O research/marcus-de-oliveira/zheng-t4-chart-sign/check.py
python3 -m unittest discover -s research/marcus-de-oliveira/zheng-t4-chart-sign -p 'test_*.py'
```

The replay checks identities in rational-function fields, all 24 raw products,
the six contacts, inverse maps, Hermitian coefficients, determinants,
evaluations, slack and arrowhead Schur complement. `replay.json` records the
reconstructed matrices. No author data, optimization, random trials or floating
point is used. The six-gap implication and eigenvalue argument are written
proofs, not Lean-formalized theorems. Tests deliberately reject altered
identities and non-Hermitian or higher-degree forms, including under `-O`.

A separate exploratory 84-point rational fixture search did not find a
seven-gap point with `d>0`; it supports no claim of emptiness and is not part of
the proof. The current second-pass audit gives numerical examples of the thin
supported `d>0` branch; an exact representative would be a useful follow-up.
