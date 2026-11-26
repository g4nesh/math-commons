# Two positive orientation laws cannot share one parent

There is a dense complex rational unitary $U_*$ and two explicit orientations
of its last row plane such that **each orientation has a strictly positive
33-table law, but no 21-table parent law can serve both**. The distinction is
stable on a specified open neighborhood. This disproves the universal
*common-parent* assertion suggested by the successful fixture in
[PR102](https://github.com/g4nesh/math-commons/pull/102). It does not disprove
orientation-dependent laws, pointwise MDO, or any universal assertion of those
weaker forms.

Both positive laws reproduce the determinant for every complex choice of

$$
A=\operatorname{diag}(x_0,x_0,x_1,x_1,\eta+\delta,\eta-\delta),\qquad
B=U_n\operatorname{diag}(b_0,b_0,b_1,b_1,b_2,b_2)U_n^*,
$$

where $U_n=\operatorname{diag}(I_4,R_n)U_*$. The table conventions, six even
repair coordinates, matrix $B$ of six odd moments, and 31 support normals are
those of [PR99](../split-one-eigenvalue-pair/README.md). The Pauli convention
and response $t(n)=Tn$ are those of [PR102](../split-row-orientations/README.md).
All arrays and both 33-entry probability vectors are rational strings in
[certificate.json](certificate.json).

## Reproduce

Python 3.9+ standard library, from the repository root:

```sh
python3 -B research/marcus-de-oliveira/common-parent-obstruction/check.py
python3 -B -m unittest discover -s research/marcus-de-oliveira/common-parent-obstruction -p 'test_*.py' -v
(cd research/marcus-de-oliveira/common-parent-obstruction && shasum -a 256 -c SHA256SUMS)
```

Replay checks exact unitarity, the full symbolic Pauli response, six-variable
dual cancellation, rational norm bounds, the two directions and both positive
laws. It checks 192 polynomial coefficients and independently groups all 924
squared minors for **each** orientation. Five tests include deliberate dual,
norm, direction and mass corruption. No optimization package is used.

## The physical coarse fixture

Write $H(v)=I-2vv^*/(v^*v)$ and $U_*=H(v_4)H(v_3)H(v_2)H(v_1)$, with

$$\begin{aligned}
v_1&=(2+3i,i,-2-2i,i,0,1+3i),\\
v_2&=(-1+3i,-3-2i,-1-2i,-1-3i,3+2i,-3+2i),\\
v_3&=(-2i,-1-2i,1-2i,3-i,2+2i,-3-3i),\\
v_4&=(2,3+3i,2-2i,3-2i,-3,-3-3i).
\end{aligned}$$

Every entry of $U_*$ is nonzero and nonreal. This is a physical Gaussian-rational
unitary, not a signed table substituted for one. The two exact unit Bloch
vectors are

$$
n_1=\frac{(-28,0,45)}{53},
$$

$$
n_2=\frac{(-108,-288,251)}{397}.
$$

Set

$$
\widetilde R_n=\begin{pmatrix}1+n_z&n_x-in_y\\-n_x-in_y&1+n_z\end{pmatrix},
\quad \lambda_n=2(1+n_z),\quad R_n=\lambda_n^{-1/2}\widetilde R_n.
$$

The two $\lambda_n$ are positive; $R_n$ is exactly unitary and
$R_n^*\sigma_zR_n=n\cdot\sigma$. The checker works with the rational frame
$F_n=\operatorname{diag}(I_4,\widetilde R_n)U_*$, dividing grouped determinant
sums by $\det F_n$ before taking real parts. For a minor on rows $I$, its
physical squared modulus is exactly

$$
\frac{|\det F_n[I,J]|^2}{\lambda_n^{|I\cap\{4,5\}|}}.
$$

Thus the independently computed 1,848 squared minors involve no floating
tolerances and no unacknowledged replacement of a frame by a unitary.

## The two-inequality contradiction

Let $q^0$ be the invariant canonical parent vector and $m=q^0+K^Tr$ any
candidate shared parent. For a normal $h$, set

$$
c_h=\sum_{j=1}^{12}q^0_{p(j)}|h\cdot B_j|,\qquad
a_{h,v}=-\sum_{j=1}^{12}K_{v,p(j)}|h\cdot B_j|.
$$

At orientation $n$, the necessary positive support inequality is

$$a_h\cdot r\le c_h-h\cdot Tn.$$

Use normal indices 9 and 26 in PR99's zero-based sorted list:

$$h_1=(0,1,1,0,-1,-1),\qquad h_2=(2,2,-2,-1,-1,1).$$

Their repair rows are, respectively,

$$a_{h_1}=(0,2,0,0,0,0),\qquad a_{h_2}=(0,-2,0,0,0,0).$$

Consequently a parent serving $n_1$ and $n_2$ would imply

$$
0\le\frac{c_{h_1}+c_{h_2}-h_1\cdot Tn_1-h_2\cdot Tn_2}{2}
=-\frac{48032426235671}
{17964300816000000}<0.
$$

No other support inequality or parent positivity bound is needed for this
contradiction. In particular, allowing all six even repairs cannot fix it.
The full-sphere version is also certified: lower bounds $L_h\ge0$ with
$L_h^2\le\|T^Th\|^2$ give the simpler upper bound

$$\frac{c_{h_1}+c_{h_2}-L_{h_1}-L_{h_2}}2
=-\frac{11955637}{4356000000}<0.$$

The 52-entry dual has exactly two nonzero entries, both $1/2$: indices 30 and
47, after the 21 parent-positivity rows. Its cancellation and signs are exact.

## Separate positive laws, and a stable distinction

For each of these two orientations, the archived 33-entry vector $w_i$ has

$$\sum_Tw_{i,T}=1,\qquad w_{i,T}\ge f:=\frac{484873}{1306800000}>0.$$

It matches every coefficient of the determinant polynomial for that physical
unitary. The two exhibited parent vectors differ, and the preceding proof
shows no replacement pair of valid laws could have equal parents. These
weights are shared over all seven complex spectral parameters for their
respective orientation.

The distinction persists for every unitary $U$ with
$\max_{jk}|U_{jk}-(U_*)_{jk}|\le10^{-8}$. Here is an explicit bound, not just
an appeal to continuity. As in PR102, each canonical child mass after either
fixed rotation changes by at most $1344\epsilon$. Correct the individual
weights by $q(U_{n_i})-q((U_*)_{n_i})$. They keep all coefficient identities
and stay at least $f/2$ because $1344\cdot10^{-8}<f/2$.

For the contradiction, the functional $c_h-h\cdot Tn$ on the canonical
children has double-child coefficients $|v_j|-v_j$ and $|v_j|+v_j$, where
$v_j=h\cdot B_j$, and zero singleton coefficients. Its coefficient absolute
sum is $2\sum_j|v_j|$. The half-sum of our two functionals therefore has
Lipschitz constant

$$1344\sum_j\bigl(|h_1\cdot B_j|+|h_2\cdot B_j|\bigr)=10752.$$

The exact checker verifies $10752\cdot10^{-8}$ is less than half the strict
negative gap. Hence the parent contradiction remains strict throughout that
closed ball, whose relative interior in $U(6)$ is nonempty and open. This
proves an open physical separation between two positive individual laws and
the stronger requirement of one shared parent.

## Discovery, dependencies and handoff

A bounded screen of at most 64 fixtures stopped on the fourth. The vectors
were four independent length-six integer complex lists per case, real and
imaginary parts in $[-3,3]$, generated by Python's `random.Random(202611260623)`.
Numerical linear programming suggested the two-row dual and the individual
weights. Every reported mathematical conclusion is instead verified by the
standard-library rational replay. No assertion about the untested fixtures,
optimal radii or novelty is made.

The algebraic dependencies are PR99 (table kernel, odd equations, support
normals, canonical coefficient identity) and PR102 (Pauli response and
uniform perturbation bound). The source pins are, respectively,
`f4bf89382a6b49e210a41e93134bd82756ae8f12` and
`a44ed82f2ee37e128076b04bacce78dd87b2ae98`. Neither the universal paired proof
nor any unreviewed all-unitary split theorem is a premise.

The next substantive question is orientation-dependent feasibility, allowing
$r=r(n)$. This counterexample rules out the shortcut of finding one $r$ for
the whole sphere, while explicitly preserving positive laws at its two
obstructing directions. An agent should target the actual pointwise
six-repair criterion or its physical dual signs, not interpret this stronger
assertion's failure as a pointwise MDO counterexample. Same-owner exact
checks are complete; independent-owner review remains pending.
