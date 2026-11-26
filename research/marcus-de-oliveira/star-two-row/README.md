# Exact two-row compression of the unitary-minor inequality

This is a partial analytic reduction for the four-dimensional inequality ★,
not a proof of ★ on all of U(4), and not a resolution of the Marcus–de Oliveira
conjecture. It gives an exact equivalence with a quadratic and a cubic, proves
★ on the entire stratum U14=0, and reduces the entire stratum U11=0 to an
explicit real-variable inequality. A complete exact unitary refutes the
shortcut of asking all coefficients of the cubic to be nonnegative.

The compression was suggested by the coordinator. The derivation and scoped
proofs here are by ganesh-reductions. Auxiliary same-owner algebra checks are
not independent mathematical review. No novelty claim is made.

## 1. Notation and quantifiers

For U in U(4), let B_ij=|U_ij|² and D[S,T]=|det U[S,T]|². Indices in this
document are one-based. Set

    A = D[12,12],              C = D[13,14],
    g0 = B11+B33+D[14,34]−D[12,14],
    f_t/t = g0+(1+1/t)A+tC             (t>0).

The inequality ★ is g0+A+2√(AC)≥0. It is equivalent to f_t≥0 for every t>0:
the infimum of A/t+tC is 2√(AC), including zero cases by a limit. When A,C>0,
the minimizing t is √(A/C). Thus a result for t≥1/5 covers A≥C/25, not every
small A. The earlier [small-minor boundary proof](../star-small-ratio/README.md)
is separate from this reduction.

Multiplying columns by phases preserves every quantity above, so take the
first row x=(u,v,w,z) real and nonnegative, u²+v²+w²+z²=1. Let a be the second
row, a complex unit vector orthogonal to x, and write

    m_ij=x_i a_j−x_j a_i,        A=|m12|².

The remaining rows b,c range over every orthonormal basis of the complex
two-plane P=span{x,a}⊥. We use P also for its orthogonal projection. No
restriction is placed on these completions.

## 2. Eliminate both remaining rows exactly

Define real vectors E=e3, W=u e4−z e1, V=w e4−z e3 and the Hermitian operator
on the two-plane

    H = P(EE*+t WW*−VV*)P restricted to range(P),
    K = u²+||PV||²−|m14|²+(1+1/t)A.

Then f_t/t=K+⟨b,Hb⟩. Indeed, C=|⟨b,W⟩|², B33=|⟨b,E⟩|² and
D[14,34]=|⟨c,V⟩|²=||PV||²−|⟨b,V⟩|². The row-inner-product convention or its
complex conjugate gives the same real quadratic value. Every unit b in the
plane admits a unit c completing the basis. Therefore, for fixed x,a,t,

    f_t≥0 for every row completion ⇔ K I+H is positive semidefinite.

For a two-by-two Hermitian matrix this is equivalent to nonnegative trace
and determinant: its eigenvalues are real, their product is nonnegative,
and a nonnegative sum rules out two negative eigenvalues. Zero and repeated
eigenvalues are included.

Everything can be written using x and a alone. Put

    k = A+u²+w²+z²−|m34|²−|m14|²,
    h = 1−w²−|a3|²−w²−z²+|m34|²,
    d = u²+z²−|m14|²,
    e = |(1−w²)a2+vw a3|²−z²|a2|².

Then K=A/t+k, tr H=h+td and det H=te−w²A. For the trace, use
||PE||²=1−w²−|a3|², ||PW||²=u²+z²−|m14|² and
||PV||²=w²+z²−|m34|². For the determinant, expand the three rank-one terms
in any orthonormal basis of the plane:

    det H = t||PE∧PW||²−||PE∧PV||²−t||PW∧PV||².

This identity follows directly from the two-by-two determinant: each pair
of rank-one terms contributes the product of its real weights times its
squared two-by-two determinant. Since x,a are orthonormal, the projected
wedge norms are the absolute values of the four-by-four determinants with
rows x,a and the indicated two real vectors. Their complex conjugates have
the same absolute values. Direct expansion and x·a=0 give, up to signs,

    det[x;a;E;W] = (1−w²)a2+vw a3,
    det[x;a;E;V] = w m12,
    det[x;a;W;V] = z a2.

Consequently the exact necessary and sufficient conditions are

    T(t) = t tr(KI+H) = 2A+(2k+h)t+d t² ≥ 0,

    F(t) = t² det(KI+H)
         = A²+A(2k+h)t+[k²+kh+A(d−w²)]t²+(kd+e)t³ ≥ 0.

This is the complete two-row reduction, not a positivity certificate for
the displayed polynomial. In particular, its coefficients need not all
be nonnegative.

## 3. A complete zero-entry stratum: U14=0

Here z=0. The original definition immediately becomes

    g0 = u²(1−|a4|²)+|b3|²+w²|c4|² ≥ 0.

Every row has norm one, so |a4|²≤1. Thus ★ and f_t≥0 for every t>0 hold
for every complex unitary with U14=0, including all further zero entries
and every row completion. No division or nondegeneracy assumption is used.

## 4. The whole U11=0 stratum reduces to a real problem

Set u=0 and q²=w²+z². If q=0, then z=0 and section 3 applies. Suppose q>0.
An orthonormal real basis of x⊥ is

    F0=e1,   G=(0,0,z,−w)/q,   J=(0,q,−vw/q,−vz/q).

Write a=αF0+βG+γJ with |α|²+|β|²+|γ|²=1. Put

    s=|α|²,    r=|γ|²,    R=|vwβ+zγ|²/q².

Substituting a2=qγ, a3=(zβ−vwγ)/q, m34=−qβ and m14=−zα into section 2 gives

    A=v²s,
    k=(1−z²)s+q²r,
    k+h=2A+R,
    d=z²(1−s),
    e=z²(R−q²r).

For example the e identity uses
(1−w²)a2+vw a3=(z/q)(vwβ+zγ). For k+h, expansion gives
v²s+1−w²−|a3|²−z²s; expanding |a3|² and using
s+|β|²+|γ|²=1 yields 2v²s+R. These are identities, not inequalities.

The trace condition is now automatic:

    T(t)=2A+(k+2A+R)t+dt² ≥ 0.

All terms are nonnegative. The determinant polynomial is

    F(t)=A²+A(k+2A+R)t
         +[k(2A+R)+A(d−w²)]t²
         +[kd+z²(R−q²r)]t³.

For fixed v,w,z,s,r,t, this expression is increasing in R, with slope

    A t+k t²+z²t³ ≥ 0.

The exact smallest permitted R is

    R_min = (vw√(1−s−r)−z√r)²/q²,

attained by choosing β=√(1−s−r), γ=−√r. This follows from the reverse
triangle inequality, and works when either magnitude is zero. Choose
α=√s as well; its phase does not enter any of these coefficients. Thus
the same real second row simultaneously minimizes F(t) for every t>0.

It follows rigorously that ★ on all complex unitaries with U11=0 is
equivalent to ★ on real orthogonal matrices with U11=0. For the nontrivial
direction, assume the real statement. For each complex x,a, take the real
row a_minus just constructed. The plane orthogonal to real x,a_minus has
a real orthonormal basis; its compressed pencil is real symmetric. Its
minimum Rayleigh quotient has a real unit eigenvector b and a real unit
completion c. The assumed real statement therefore makes its determinant
nonnegative for every t: the minimum Rayleigh quotient of K I+H is
nonnegative, so the matrix is positive semidefinite. Monotonicity in R and
the automatic nonnegative trace imply positivity for the original complex
x,a and every completion.
The reverse implication is inclusion of real matrices in complex matrices.

This is a reduction on one zero-entry stratum, not a reduction of all U(4)
to O(4), and not a proof of the remaining real statement.

### Explicit remaining low-dimensional inequality

For every v,w,z≥0 with v²+w²+z²=1 and q²=w²+z²>0, every s,r≥0 with s+r≤1,
and every t>0, prove F(t)≥0 with the formulas immediately above and
R=R_min. There are two independent first-row parameters, two magnitude
parameters and t; no complex phases or row-completion variables remain.
The zero-dimensional q=0 boundary is already proved in section 3.

There is also a radical-free form. Let B,G≥0 with B²+G²≤1, set

    Q=w²+z²,   S=1−B²−G²,   A=v²S,
    K=(1−z²)S+QG²,   D=z²(1−S),   L=(vwB−zG)².

The remaining statement, for Q>0, is precisely nonnegativity of the
following polynomial for every admissible v,w,z,B,G,t:

    Q A² + A[Q(K+2A)+L]t
      +[K(2QA+L)+QA(D−w²)]t²
      +[QKD+z²(L−Q²G²)]t³ ≥ 0.

This is QF(t), with s=S, r=G² and R=L/Q; the positive multiplier preserves
the sign. No sampling assumption enters this equivalence.

## 5. Exact obstruction and reproducible checks

`negative-c2.json` preserves a complete exact Gaussian-rational unitary
and its four cubic coefficients. Its t² coefficient is

    −6967634414959095922166053753681234854015302933504
    /30795100483993887135129899427539005627612941537375625.

`u11-negative-c2.json` has U11=0 and coefficient

    −1426840037349781197994246969344
    /57449173098767439652925611224715625.

Both are strictly negative exact rational numbers. The replay script checks
unitarity and recomputes every coefficient from the saved matrices. This
refutes coefficientwise positivity, including on U11=0. It is not a negative
value of F(t), a counterexample to ★, or evidence of a universal sign for
any other coefficient.

From the repository root, with Python 3 and no external packages:

```sh
python3 -B research/marcus-de-oliveira/star-two-row/verify.py
python3 -B research/marcus-de-oliveira/star-two-row/scan.py
```

`verify.py` checks the compression directly in the actual last-two-row basis
at t=1/100,1/5,1,5, recomputes both saved obstructions, and checks the
zero-stratum identities and phase reduction. `scan.py` exactly reproduces
the fixed-seed bounded discovery runs: 1500 unrestricted-first-row inputs
and 1000 U11=0 inputs. These parametrized finite samples are not exhaustive
over U(4). Reports distinguish exact verification from exploratory counts.
`--write-reports` is optional and rewrites only the local report/fixture files.
Neither script uses floating-point arithmetic to decide an equality or sign.

The arithmetic/frame helper is imported read-only from
`../star-small-ratio/verify.py`, pinned at source commit
`3a9f06ba02a650597879f17a47da7bccb966481b`; its SHA-256 is checked before use.
The full matrix is stored as pairs of rational strings, never rounded floats.
The analytic proofs above carry the universal quantifiers; finite replay
supports the algebra and the explicit coefficient obstruction only.

The saved replay used Python 3.9.6: 384 general cases, 384 U14=0 cases,
384 U11=0 cases, 384 opposite-real phase comparisons, 24 permutation
degeneracies, and both negative-coefficient fixtures. Verification used
about 11 process CPU seconds; the two discovery replays used about 5 seconds
combined. The respective scans have 4/1500 and 2/1000 negative middle
coefficients. These counts describe this exact generator, not a frequency
under Haar measure or an exhaustive search.

## 6. Remaining work

The global quadratic/cubic inequalities and the final real U11=0 polynomial
are open in this artifact. A successful route must control cancellation in
the cubic rather than assert coefficientwise positivity. The two-row
compression is exact for every t>0 and includes degeneracies, so a negative
eigenvalue would provide a concrete row completion; none is asserted here.
This work contains no Lean formalization or independent-owner review.
