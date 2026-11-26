# The 45 Cauchy–Binet inequalities do not imply pointwise spectral support

The proposed universal cone inclusion is **false**, even for explicit Gaussian-integer spectra. This is an obstruction to using only the 45 inequalities, not a counterexample to Marcus–de Oliveira: the separating signed table below is not asserted to arise from a unitary.

The starting reduction is [PR81](https://github.com/g4nesh/math-commons/pull/81), with [immutable source](https://github.com/g4nesh/math-commons/blob/272f7ced73c3b880cff0f44adfe25f6f8cceca7f/research/marcus-de-oliveira/double-paired-contingencies/README.md). That work identifies 21 unordered pairs of $S_3$ permutations, a six-dimensional polynomial repair kernel, and 189 dual rays in seven symmetry orbits. The first three orbits contain 45 rays supplied by grouped Cauchy–Binet positivity. This note gives one exact normalized dual vertex passing those 45 inequalities but failing the requested **pointwise** spectral support test.

## Conventions and the cone

Use zero-based permutations in the order

$$\pi_0=012,\quad\pi_1=120,\quad\pi_2=201,\quad
\pi_3=021,\quad\pi_4=102,\quad\pi_5=210.$$

The first three are even and the last three odd. Put $p_i=\prod_{r=0}^2(a_r+b_{\pi_i(r)})$. Coordinates are the 21 pairs $(i,j)$ with $0\le i\le j<6$, in lexicographic order; no off-diagonal factor is implicit in the pairing $q\cdot G$.

For an even–odd edge $e=\{a,b\}$, write $g_e$ for its six-dimensional incidence vector. For the 45 unordered pairs of the nine edges, define

$$R_{ef}(i,j)=g_e(i)g_f(j)+g_f(i)g_e(j),\qquad K=\operatorname{cone}\{R_{ef}:e\le f\}.$$

The diagonal choice $e=f$ therefore includes a factor of two. Rescaling these generators gives exactly the first 45 primitive rays of PR81. The candidate claim was that, for every spectrum and complex direction $\theta$,

$$G_{ij}=\Re(\theta p_ip_j)-\min_{k\le l}\Re(\theta p_kp_l)\quad\text{belongs to }K.$$

## Exact obstruction

Take

$$a=(6-i,\;1+6i,\;-6-8i),\qquad b=(-6-8i,\;-1+2i,\;3-i),\qquad\theta=-1+2i.$$

For the two-sided $(2,2,2)$ problem, each entry of each triple is repeated twice. Their six pairing products are

$$p=(-216-648i,\;284-588i,\;295+350i,\;-531-18i,\;-66+252i,\;960-1120i).$$

Define the 21-coordinate signed table, in the specified lexicographic order, by

$$4q=(2,2,2,-1,0,-2,\;0,1,0,1,0,\;1,-1,0,0,\;0,-1,1,\;-1,0,\;0).$$

The checker computes all of the following with integers and fractions:

* $\sum q_{ij}=1$.
* Every one of the 45 numbers $q\cdot R_{ef}$ is $0$, $1/4$, or $1/2$.
* Every one of the 141 grouped monomial coefficient rows also pairs nonnegatively with $q$.
* The minimum vertex score is $-850560$, attained only at $(i,j)=(4,5)$.
* The abstract signed value is $\sum q_{ij}p_ip_j=120416+375488i$, whose directional score is $-871392$.

Consequently

$$q\cdot G=-871392-(-850560)=-20832<0.$$

If $G$ were a nonnegative combination of the 45 generators, pairing with $q$ would be nonnegative. This contradiction proves $G\notin K$ exactly. All the pairing values come from actual spectra; the only relaxed object is the signed table $q$.

## It is a vertex of the normalized dual quotient

The 45 generators have exact rank 15 in $\mathbb R^{21}$. Their annihilator is the six-dimensional polynomial repair kernel from PR81. A dual table is considered modulo this annihilator. Normalization $\sum q=1$ is well defined on the quotient: the all-ones row is a polynomial coefficient row and annihilates that kernel.

There are exactly 20 active generator inequalities at this $q$, and their rows have rank 14. Adding the normalization row raises the rank to 15. Thus the active constraints and normalization determine one point of the 15-dimensional dual quotient. Any convex decomposition into feasible normalized tables forces both summands to satisfy every active equality, hence to be the same quotient point. This proves the vertex claim without enumerating all vertices.

For completeness, the normalization slice is bounded in the quotient. If $S=\sum_e g_e=3\mathbf1$, then

$$\mathbf1\mathbf1^T=\frac19\sum_{e<f}R_{ef}+\frac1{18}\sum_eR_{ee}.$$

All coefficients are positive, so normalization is strictly positive on every nonzero ray of the pointed dual quotient.

## Precisely which additional inequalities repair this example

The table violates seven of PR81's 144 additional rays. In the following display, $[ij]$ denotes the coordinate row selecting $(i,j)$. Orbit numbers refer to the representatives in PR81, not to a new classification.

| Extra ray | PR81 orbit | Pairing with $q$ |
| --- | --- | --- |
| $2[02]+2[05]+[23]+[24]+[34]+2[55]$ | 5 | $-1/2$ |
| $[02]+[05]+[24]+[44]+[55]$ | 4 | $-1/4$ |
| $[02]+[05]+[23]+[33]+[55]$ | 4 | $-1/4$ |
| $[02]+[05]+[11]+[14]+[23]+[34]+[55]$ | 6 | $-1/4$ |
| $[02]+[05]+[11]+[13]+[24]+[34]+[55]$ | 6 | $-1/4$ |
| $[01]+[02]+2[05]+[14]+[23]+[34]+2[55]$ | 7 | $-1/4$ |
| $[01]+[02]+2[05]+[13]+[24]+[34]+2[55]$ | 7 | $-1/4$ |

The checker expands the seven published symmetry orbits and verifies that this list is exhaustive. The violation counts in orbits 1–7 are $(0,0,0,2,1,2,2)$.

More constructively, **each single extra ray in this table separately suffices to repair this particular spectral gap**. For every listed $r$, [certificate.json](certificate.json) supplies an exact identity

$$G=\sum_{e\le f}\lambda_{ef}R_{ef}+\mu r,\qquad\lambda_{ef}\ge0,\quad\mu>0,$$

with 14 positive terms total. The standalone checker verifies all seven identities coordinate by coordinate. This is stronger than merely observing that $q$ violates the extra inequalities. It does not prove any of those inequalities for physical unitaries or prove that one orbit suffices for every spectrum.

## Reproduction, evidence level and next step

Python 3.9+ standard library, from the repository root:

```sh
python3 -B research/marcus-de-oliveira/paired-pointwise-cone/check.py
python3 -B -m unittest discover -s research/marcus-de-oliveira/paired-pointwise-cone -p 'test_*.py'
```

The checker reconstructs the 45 generators, all 141 coefficient rows, the 189 published orbit rays, exact vertex ranks, spectral products, strict separation and seven positive repairs. It compares [report.json](report.json). It has no LP solver or dependency on another research directory. Discovery used a bounded SciPy LP probe; the obstruction appeared after 1,795 actual spectrum samples in that run. Discovery samples are not the proof. The exact certificate is.

This result closes the proposed 45-inequality sufficiency route negatively. It does not claim a physical MDO violation, an independent-owner review or a Lean proof. No unitary realizing $q$ is supplied. It does not duplicate the separate contact-reduction theorem, and the example initially has one exposed contact. A concrete next investigation is whether physically valid constraints from any one of the four missing PR81 orbits suffice for all spectral support directions; all four can repair this particular example, but universal coverage remains open.
