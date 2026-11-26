# Which contact patterns at n = 5 are automatic?

Task `mdo-n5-contact-strata` · agent `ishaan-moc` · owner IshaanAyaan · 26 November 2026.

**Deviation from the task text.** The task asks for sampling by the contact deformation. On the coordinator's instruction this lane works exact-first: it reduces the question to the facets of one polytope and certifies every facet exactly. The deformation sampling has not been run, so nothing here decides which contact sets are *realizable* as rigid contacts.

## Result so far

1. **Reduction (proved below).** A contact set $T$ is *automatic*, meaning that every generator of its gap cone is a squared-minor vector, if and only if no facet of the minor polytope $Q_n$ that contains $T$ is a non-minor facet.
2. **Check at n = 4 (EXACT).** With this criterion, own code reproduces the dimension-four classification exactly:
   - $Q_4$ has 82 facets, 34 squared-minor and 48 others;
   - the 6-element sets fall into 214 orbits, 124 automatic and 90 not.
3. **n = 5.**
   - $Q_5$ has 120 vertices and dimension 41. Seven orbits of non-minor facets are certified exactly (table below). Whether this list is complete is **not** known.
   - **Every tested 8-set is non-automatic (EXACT).** An exact cutting-plane census ([`census.py`](census.py)) repeatedly asks an integer program for an 8-set that meets the support of every known non-minor facet. These are exactly the 8-sets that could still be automatic. It then decides each candidate in exact rational arithmetic. All 9 candidates produced so far are non-automatic:
     - 2 by an exact linear program that finds a gap vector outside the span of the squared-minor vectors;
     - 7 by exact enumeration of all extreme rays of their gap cone, which have between 7 and 419 non-minor rays.
     Five of these decisions supplied the new orbits C–G below; the remaining non-minor rays are recorded in `census_n5.json`.
   - **Conjecture.** No 8-element set of permutations is automatic at $n=5$: every one lies on a non-minor facet of $Q_5$. Two pieces of evidence:
     - all 200,000 random 8-sets lie on a facet of orbit A alone;
     - after the seven orbits, HiGHS finds no 8-set meeting all 57,780 non-minor supports within 600 s, nor within 1,500 s when the 8-set is required to contain the identity (34,122 supports avoid it).
     A solver time-out is not a proof.

| Non-minor facet orbit | Vertices on the facet | Values of the gap vector off the facet | Orbit size (group of order 28,800) | Found by |
| --- | --- | --- | --- | --- |
| A | 84 | 1 on 36 | 900 | random LP search |
| B | 74 | 1 on 44, 2 on 2 | 14,400 | random LP search |
| C | 43 | 1 on 45, 2 on 25, 3 on 5, 5 on 2 | 2,880 | census |
| D | 74 | 1 on 44, 2 on 2 | 14,400 | census |
| E | 76 | 1 on 40, 2 on 4 | 3,600 | census, ray enumeration |
| F | 59 | 1 on 51, 2 on 9, 3 on 1 | 7,200 | census, ray enumeration |
| G | 74 | 1 on 44, 2 on 2 | 14,400 | census, ray enumeration |

B, D and G have the same values but are distinct orbits (different canonical vertex sets).

Squared-minor facets are $x_{K,J}\geq0$:
- $|K|=1$: 96 vertices;
- $|K|=2$: 108 vertices.

The 82 facets of $Q_4$ and the seven $n=5$ orbits above are all certified the same way: the vertices on each facet have affine rank $\dim Q_n$, and the affine function vanishing on them is strictly one-signed on the remaining vertices, in exact rational arithmetic (python-flint).

**What this means for the program.** At $n=4$, 124 of 214 patterns were automatic. At $n=5$, the squared-minor cone settles no 8-contact pattern tested so far, and conjecturally none at all. Nearly every rigid configuration will need a case argument, as the 90 non-automatic cases did at $n=4$. So a classification "up to minors" will not shrink the $n=5$ problem the way it did at $n=4$.

## Structure of the non-minor facets

Write $x_{K\to J}(\sigma)=1[\sigma(K)=J]$ (with $x_{K\to J}=x_{K^c\to J^c}$). Minimum-$\ell_1$ representations of the seven gap vectors in these coordinates, with integer coefficients verified exactly on all 120 vertices, are:

| Orbit | Facet inequality (one representative) |
| --- | --- |
| A | $x_{12\to12}\le x_{0\to0}+x_{03\to34}+x_{04\to34}$ |
| B | $x_{12\to24}\le x_{1\to2}+x_{02\to23}+x_{03\to01}+x_{13\to34}$ |
| D | $x_{23\to23}\le x_{03\to34}+x_{04\to04}+x_{12\to24}+x_{14\to04}+x_{24\to13}$ |
| E | $x_{03\to04}\le x_{2\to2}+x_{12\to13}+x_{12\to23}+x_{14\to23}$ |
| F | $x_{04\to23}\le x_{3\to0}+x_{12\to01}+x_{14\to03}+x_{14\to34}+x_{24\to02}+x_{24\to24}$ |
| G | $x_{04\to03}\le x_{02\to04}+x_{03\to23}+x_{12\to12}+x_{12\to24}+x_{14\to34}$ |
| C | no integer representation with few terms; the minimum-$\ell_1$ form has 25 fractional terms |

Six of the seven are **covering inequalities**: a single squared-minor event is bounded by a sum of others that together cover every permutation on which it occurs. Validity is a finite counting statement.

Orbit A is a block inequality. Take a row 3-set $S=\{a,b,c\}$ and a column 3-set $S'=\{a',b',c'\}$; $x_{\bar S\to\bar S'}=x_{S\to S'}$. Each of the six bijections $S\to S'$ makes exactly one of $x_{a\to a'}$, $x_{ab\to b'c'}$, $x_{ac\to b'c'}$ equal to 1, so equality holds on all 12 permutations with $\sigma(S)=S'$.

These facets are inequalities on the vertices of $Q_5$ only. They are not claimed for the squared minors of a unitary matrix.

## The reduction

Notation:
- $v_\sigma\in\{0,1\}^{\mathcal P}$ has coordinates indexed by the pairs $(K,J)$ with $1\leq|K|=|J|\leq n-1$, identified with $(K^c,J^c)$. The entry is $1$ iff $\sigma(K)=J$.
- $Q_n=\mathrm{conv}\{v_\sigma\}$.
- $L\subset\mathbb R^{S_n}$ is the space of vectors $(\ell(v_\sigma))_\sigma$ for affine functions $\ell$. It is the span of the squared-minor vectors $n_{K,J}=(1[\sigma(K)=J])_\sigma$ and of $\mathbf 1$; at $n=5$ it has dimension 42.
- $C=L\cap\mathbb R^{S_n}_{\geq0}$ is the gap cone.

For a contact set $T$, the cone in the task is $F_T=\{g\in C:g_T=0\}$.

**Lemma.** $F_T$ is generated by the vectors $(\ell(v_\sigma))_\sigma$ for the facet inequalities $\ell\geq0$ of $Q_n$ whose facet contains every $v_\sigma$ with $\sigma\in T$. Hence every generator of $F_T$ is a squared-minor vector iff each such facet is a coordinate facet $x_{K,J}\geq0$.

*Proof.*
1. The map $\ell\mapsto(\ell(v_\sigma))_\sigma$ identifies affine functions modulo those vanishing on $\mathrm{aff}\,Q_n$ with $L$. It carries the cone of affine functions nonnegative on $Q_n$ onto $C$.
2. That cone is pointed modulo the equations, and its extreme rays are the facet inequalities, by the Minkowski–Weyl theorem.
3. Because $g\geq0$, $F_T$ is the face of $C$ cut out by $g_\sigma=0$ for $\sigma\in T$. The extreme rays of a face are the extreme rays of the cone that lie in it. Hence $F_T$ is generated by the facet inequalities vanishing on $T$.
4. A facet inequality gives a squared-minor vector exactly when it is a coordinate inequality $x_{K,J}\geq0$. A non-coordinate facet's vector cannot be a positive multiple of a squared-minor vector, since distinct facets have distinct vertex sets. ∎

**Corollary.** $T$ is automatic iff $T$ meets the support of the gap vector of every non-minor facet. A superset of an automatic set is automatic.

## Files

| File | What it does |
| --- | --- |
| [`minor_polytope.py`](minor_polytope.py) | Builds $Q_n$, enumerates the facets of $Q_4$ exactly (cddlib, GMP rationals) and the orbits of 6-sets. `python3 -I -B minor_polytope.py 4` |
| [`facet_search.py`](facet_search.py) | Proposes facets of $Q_n$ by random linear programs with forced zeros, certifies each exactly and reduces to orbits. `python3 -I -B facet_search.py 5 3000 2` (about 15 minutes) |
| [`automatic_certificate.py`](automatic_certificate.py) | Exact test of one 8-set: dual certificates, dimension count and private permutations. |
| [`census.py`](census.py) | The cutting-plane census. `python3 -I -B census.py 15 600` extends `facets_n5.json` and `census_n5.json`. |
| [`census_fixed_identity.py`](census_fixed_identity.py) | Integer program for an 8-set containing the identity that meets every known non-minor support. |
| [`resolve.py`](resolve.py) | Exact extreme-ray enumeration of the gap cone for candidates whose squared-minor vectors lack private permutations. |
| [`exact_search.py`](exact_search.py) | Exact depth-first search for an 8-set containing the identity that meets every known non-minor support (does not finish; see Open). |
| [`facets_n5.json`](facets_n5.json) | All certified facet orbits at $n=5$ (2 squared-minor, 7 others), with one vertex set each. |
| [`census_n5.json`](census_n5.json) | Every census candidate and its exact verdict. |

Requirements: Python 3, numpy, scipy (HiGHS), python-flint, and pycddlib built against cddlib. Direct exact enumeration of $Q_5$ by cddlib did not finish within an hour and was stopped.

## Open

- **The conjecture.** Prove that no 8-set is automatic. A direct exact search ([`exact_search.py`](exact_search.py)) over 8-sets containing the identity, with a disjoint-packing bound, was run twice for 20 minutes (44 and 20 million nodes) and closed none of its 36 top-level branches, so brute force on a laptop is not a route. One route is an exact certificate that the hitting-set program is infeasible: either a branch-and-bound tree with exact bounds, or a symmetric exhaustive search over the 8-sets containing the identity. Another is a structural description of the non-minor facets of $Q_5$.
- **Completeness of the facet list of $Q_5$.** Not needed for the conjecture. Any non-automatic verdict only needs one non-minor facet containing $T$.
- **Realizability.** Which patterns occur as rigid contact sets requires the deformation of PR 71, which is not run here. If the conjecture holds, every rigid configuration at $n=5$ needs its own case argument.
