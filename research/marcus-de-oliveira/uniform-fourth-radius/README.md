# An explicit fourth-eigenvalue radius for every complex triple

Fix the archived real orthogonal four-cycle matrix `spectral.U` and
\(b=(3-2i,2,1+2i,-i)\). Put \(r=10^{-100}\). For every
\[
\max(|\Re\delta|,|\Im\delta|)\le r,\qquad x,z,y\in\mathbb C,
\]
there are **24 strictly positive real probabilities**, independent of \(x,y\), such that
\[
\det\!\left(\operatorname{diag}(x,z,y,\delta)+U D_bU^T\right)
 =\sum_p m_p\prod_{j=1}^4(a_{p(j)}+b_j),\qquad a=(x,z,y,\delta).
\]
Every probability is greater than \(1/2420000000\). The law may depend on \(z,\delta\).
This quantifies the [projective compactness theorem](../projective-fiber-stability/uniform/README.md)
at its fixed fixture. It does **not** provide a numerical neighborhood in \(U,b\),
claim an optimal radius, or resolve universal MDO. Independent-owner review is pending.

## Equations and a signed reference law

Write \(h=x-\delta\), \(k=y-\delta\), \(w=z-\delta\), and \(b'=b+\delta\mathbf1\).
For a row-to-column word \(q=(i,j,l,t)\), the four complex raw coefficient entries are
\[
C=b'_i b'_l b'_t,\quad L=(w+b'_j)b'_l b'_t,\quad
M=(w+b'_j)b'_i b'_t,\quad N=(w+b'_j)b'_t.
\]
Together with the mass row, their real and imaginary parts form a real
\(9\times24\) matrix \(A(w,\delta)\). Its right side \(d(w,\delta)\) contains
normalization and the corresponding determinant coefficients. The constant
term minus the common value at \(h=w=k=0\) has a factor \(w\); cancel that factor
before defining \(C\), including at \(w=0\). The other entries are coefficients
of \(h,k,hk\). A solution certifies all complex \(h,k\) simultaneously.
On the reciprocal chart \(v=1/w\), keep \(C\) unchanged and replace each
\((w+b'_j)\) in the other entries by \((1+vb'_j)\), scaling their target entries too.

We use the signed permutation law discussed by `ganesh-claude-minors` in the
[shared minors report, §4, canonical-formula remark](../minors/report.md#4-the-three-eigenvalue-conjecture).
For this **real** orthogonal fixture it also follows directly from
\[
\det(A+UD_bU^T)\det U=\det(AU+UD_b).
\]
Expanding the right side by columns gives the constant real signed masses
\[
s_p=\det(U)\operatorname{sgn}(p)\prod_j U_{p(j),j}.
\]
They sum to one, reconstruct every complex diagonal spectrum, and satisfy
\(\sum_p|s_p|\le24\), because each entry of an orthogonal matrix has modulus at most one.
The exact fixture value is \(10718023/9765625\), but the proof uses only 24.
Consequently \(A(w,\delta)s=d(w,\delta)\) in both charts, including the canceled
coefficient by polynomial identity. These masses are **not** asserted positive.
The analogous complex-unitary formula need not be real; no such extension is used.

The angular laws and signed masses use column-to-row words \(p\); the rank
certificates use their inverses \(q=p^{-1}\). The checker explicitly converts them.

## A common positive floor before moving the fourth entry

Replay all eight corrected angular laws from the projective theorem. In a patch
with even radial degree \(R\in\{2,4\}\), every numerator contains
\(10^{-7}(1+u^R)\), including at every Bernstein control index. Adjusted archived
coefficients remain positive and sum to one at each control. The two degree-four
patches have radial centers of modulus at most 2, while the six degree-two
patches have centers of modulus at most **10**.

If \(|u|\le1\), each degree-two radial square is at most \(11^2=121\), and
each degree-four product is at most \(3^4=81\). If \(|u|\ge1\), the same
bounds hold after division by \(|u|^R\), also for lower-degree terms. Thus the
normalizing denominator is at most \(121\max(1,|u|^R)\). Every base mass is at least
\[
\mu=\frac{10^{-7}}{121}=\frac1{1210000000}.
\]
Bernstein coefficients are nonnegative and sum to one, so angular endpoints
are included. Radial limits at infinity have the same bound. At \(w=0\),
the canceled \(C\) equation follows by continuity of a radial law. This gives
one suitable full-support law \(m^0(w)\) for every projective parameter; a
single globally continuous choice is not needed.

## A quantitative cover of the whole middle-eigenvalue plane

Use \(T=25\), \(s=1/4096\), \(\eta=s/1000\), and \(\gamma=\eta^5\).
The cover consists of four closed pole boxes
\(\max(|\Re(w+b_j)|,|\Im(w+b_j)|)\le s\), the remaining bounded region
\(|w|\le T\), and the reciprocal region \(|w|\ge T\).

### Bounded region away from the poles

Write \(w=X+iY\). The replayed quadratic Bézout multipliers express
\[
g_1=6X+5Y^4-20Y^2+3Y+12,\qquad
g_2=Y^2(Y-2)(Y-1)(Y+2)
\]
as combinations of the six archived nine-column minors \(D_n\).
The sums of absolute polynomial coefficients of both multiplier families are
less than one, checked exactly. Therefore on \(|X|,|Y|\le T\),
\(\sum_n|P_n(X,Y)|\le T^2\), and likewise for the second identity.

If \(Y\) is at least \(\eta\) from every root \(-2,0,1,2\), each of the five
factors of \(g_2\), counting multiplicity, has modulus at least \(\eta\).
Thus \(|g_2|\ge\gamma\). Otherwise choose a root \(Y_j\) within \(\eta\).
Since we are outside its pole box, \(|X-X_j|\ge s\). On the short interval
between \(Y,Y_j\), \(|Y|\le3\), and
\[
|20Y^3-40Y+3|\le663<700.
\]
Hence \(|g_1|\ge6s-700\eta\ge\gamma\). The corresponding \(X_j\) is the
root of \(g_1(X,Y_j)\), precisely the real part of \(-b_j\).
It follows that at least one of the six minors satisfies
\[
|D_n|\ge d_F:=\gamma/T^2.
\]
This is an inequality on a continuous region, not a grid approximation.

### Reciprocal region

For the archived infinity basis, direct exact replay gives
\[
D_\infty(X,Y)=25000(18X-Y+4)
 (X^2-4X+Y^2-1)(2X^2-X+2Y^2+Y-1),\qquad v=X+iY.
\]
There are only six affine rows, so its total degree is at most six. Equality
on the 28-point integer grid \(X,Y\ge0\), \(X+Y\le6\), proves this polynomial
identity by triangular Newton interpolation. The sums of absolute power
coefficients by degree are
\[
[100000,1075000,3600000,4500000,9050000,3800000,0].
\]
For \(|v|\le1/T\), both real coordinates are at most \(1/T\), giving
\[
|D_\infty(v)-100000|\le\sum_{d=1}^6 c_dT^{-d}
 =\frac{153348616}{3125}<50000.
\]
Thus \(|D_\infty|>d_I:=50000\) throughout this region, in every direction.

## Moving the fourth entry without losing positive masses

All \(|b_j|\le4\). In the finite region at \(\delta=0\), every real matrix entry
has modulus at most \(K_F=100T\). If \(|\delta|\le1\), telescoping the products
above bounds their changes by \(61|\delta|\) for \(C\),
\((9T+61)|\delta|\) for \(L,M\), and \((T+9)|\delta|\) for \(N\).
These are all at most \(100T|\delta|\). On our closed square \(|\delta|\le2r\),
so every entry of \(\Delta A\) is at most \(E_F=200Tr\).
In the reciprocal region \(|v|\le1\), the equally conservative bounds
\(K_I=1000\) and \(E_I=2000r\) follow from the same product telescoping.
The mass row does not change.

For either region select the certified nine-column block \(B\), and set
\[
R=\frac{9\,8!K^8}{d}.
\]
The adjugate formula gives \(\|B^{-1}\|_\infty\le R\), because each cofactor
has at most \(8!\) products of eight entries. The perturbation of this block
has row norm at most \(9E\). The exact checker verifies
\[
9RE<\tfrac12,\qquad 50RE<\mu/2
\tag{1}
\]
for \(r=10^{-100}\), with the finite and reciprocal constants above.
Neumann inversion consequently gives \(\|(B+\Delta B)^{-1}\|_\infty<2R\).

Hold the other 15 masses at their base values. The residual of the base law is
\[
A(w,\delta)m^0-d(w,\delta)=\Delta A\,(m^0-s).
\]
This identity uses both \(A(w,0)m^0=d(w,0)\) and the signed reference law;
it does not assume a constant row span or ignore the changing target.
Since \(\|m^0-s\|_1\le1+24=25\), its maximum entry is at most \(25E\).
Solving with the selected block changes each adjustable mass by less than
\(2R\cdot25E=50RE<\mu/2\). All masses remain greater than \(\mu/2\), and
the mass equation keeps their sum equal to one. The other eight rows give the
complete coefficient identity for all first and third eigenvalues.

### Pole boxes and boundary cases

Inside a pole box, the actual middle entry is \(z=w+\delta\), so its offset
\(z+b_j\) has each real coordinate bounded by \(s+r<1/1024\). The fourth entry
also lies in that square. The [four-collision theorem](../projective-fiber-stability/README.md)
therefore provides all 24 probabilities greater than \(1/4000>\mu/2\).
This covers the pole centers and every boundary. The regular construction
covers \(w=0\); the two regular regions overlap at \(|w|=25\). All inequalities
in (1) are strict even at the closed fourth-entry square boundary. No restriction
on singularity, repeated eigenvalues, or the values of \(x,y\) was imposed.

## Exact replay and provenance

~~~sh
python3 -B research/marcus-de-oliveira/uniform-fourth-radius/check.py --check research/marcus-de-oliveira/uniform-fourth-radius/certificate.json
python3 -B -m unittest discover -s research/marcus-de-oliveira/uniform-fourth-radius -p 'test_*.py'
~~~

The standard-library replay checks the inherited eight coefficient identities,
rank Bézout identities, four quantitative pole bounds, fixture agreement,
radial normalization bound, reciprocal determinant identity, signed masses,
and every scalar inequality above. It also constructs exact moved laws at
zero, an ordinary complex parameter, a reciprocal parameter, and all four pole
box boundaries, comparing their original dense determinants. These seven
regressions guard conventions; the continuous theorem rests on the identities
and bounds, not the samples. The command has a 60-second CPU limit.

The signed-law derivation builds on Claude's shared report; the exact signed-mass
replay was also checked by `ganesh-search`. The reciprocal determinant
factorization was supplied and independently checked by `ganesh-coordinator`.
The radial denominator bound explicitly includes the centers of size 10 in
the quadratic patches. This remains a same-owner result awaiting independent
mathematical review.
