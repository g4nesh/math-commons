# Adversarial audit: stronger statements and hidden hypotheses

Task: `mdo-adversarial` · Agent: `vn-claude` (GitHub owner `vignesh-nagarajan-vn`) · Session: November 25, 2026, 16:20–16:50 UTC.

**Evidence level:** exact computation for two finite witnesses, written review notes for the reductions, and seeded floating-point spot checks. Nothing here proves or refutes the original Marcus–de Oliveira conjecture (MOC). No independent review of this artifact has happened yet.

## 1. The failed polygon enclosure is reproduced, and it does not touch MOC

The superseded preprint DMUC p2611 v1 stated Conjecture 7: $\Delta(A,B)\subseteq\operatorname{cl}J(A,B)$. Here $\Delta$ is the set of values $\det(D_a+UD_bU^*)$, and $J$ is the bounded region inside the outer boundary of the graph on the σ-points. That graph joins $z_\sigma$ and $z_\tau$ whenever $\sigma,\tau$ differ by one transposition. Kovačec's v2 (file SHA-256 `5221b16e…fe54`, the same hash as [`../literature/sources.json`](../literature/sources.json)) withdraws the conjecture and reports 3×3 counterexamples from Goyal and from Zhang and Ma.

[`polygon_enclosure.py`](polygon_enclosure.py) checks both witnesses exactly over $\mathbb Q(i)$. It reuses only the transcribed inputs in `../certificates/fixtures/` and shares no code with `../certificates/verify.py`.

| Witness | $\det(D_a+UD_bU^*)$ | Inside the MOC hull? | Outside $\operatorname{cl}J$? |
| --- | --- | --- | --- |
| Goyal | $2/3+2i/9$ | **Yes.** All facet inequalities hold, with minimum cross-product slack $2$. | **Yes.** The ray $p+t(-1,-1)$, $t\ge0$, meets no graph edge. |
| Zhang–Ma | $(2864+1840i)/49$ | **Yes.** Minimum slack $85872/49$. | **Yes.** The ray $p+t(0,-1)$, $t\ge0$, meets no graph edge. |

**Why the ray is a certificate.** The ray test is exact. If a closed ray from $p$ meets none of the closed edge segments, then $p$ lies in the unbounded component of the complement of the graph. That puts $p$ outside $\operatorname{cl}J$. The script also checks that $p$ is on no edge. Separately, it traces the outer face of the segment arrangement as a display aid and finds winding number $0$ for both points.

**Cross-checks against the source:**

- The traced outer boundaries contain Kovačec's reported "inflex" points: $(1,0)$ for Goyal and $(7672+4424i)/155$ for Zhang–Ma.
- Both determinants equal the values printed in the paper.
- Both unitaries pass $U^*U=I$ exactly.
- Kovačec says the Zhang–Ma point "seems to lie on $[z_2,z_4]$". It does not lie on that edge or on any other.

**Conclusion.** Both witnesses refute the stronger Conjecture 7 and are consistent with MOC. Results are saved in [`polygon-enclosure-results.json`](polygon-enclosure-results.json).

```sh
# from research/marcus-de-oliveira/
python -I adversarial/polygon_enclosure.py certificates/fixtures/goyal-2026.json certificates/fixtures/zhang-ma-2026.json
```

## 2. The surviving v2 theorem carries a hypothesis nobody has established

Kovačec v2 replaces Conjecture 7 with an unnumbered theorem on p. 5: **if** $\Delta(A,B)$ is simply connected, then $\operatorname{cl}J(A,B)\subseteq\Delta(A,B)$.

- The argument is sound only under that hypothesis. Each edge lies in $\Delta$ by the generalized-Givens corollary. I spot-checked that corollary's determinant formula numerically, as described in §4.
- Simple connectedness is **not** known in general. The same paper cites Bebiano's non-simply-connected example for $n=12$, and its closing questions ask whether every $\Delta$ is simply connected.
- Any later use of $J$ as a guaranteed subset of $\Delta$ must either prove simple connectedness for the relevant family or keep the result conditional.
- The direction relevant to MOC, $\Delta\subseteq\operatorname{conv}$, gains nothing from $J$.

## 3. Audit of `../reductions/README.md`, §§1–4

I re-derived each step by hand. I found **no mathematical error**. Specifically:

- The orbit equivalence holds.
- Both directions of the support-function criterion (1) hold, including point and segment hulls.
- Every symmetry-table row holds, including the $SU(n)$ normalization by $R=\operatorname{diag}((\det U)^{-1},1,\dots)$.
- The identity $X^{-1}(X+Y)Y^{-1}=X^{-1}+Y^{-1}$ in the reciprocal transform holds.
- $C\big((I-U)(I+U)^{-1}\big)=U$.
- The continuity argument for a rational witness holds.

Two **scope cautions** should be carried into any later work that uses these reductions:

1. **Direction normalization and shifts leave structured families.** The θ=0 reduction rescales both spectra by $e^{-i\theta/n}$. It is valid for the universal problem over all spectra, as the file says. It does not hold inside a restricted class:
   - Hermitian $A$ (real spectrum $a$) acquires a nonreal spectrum.
   - Shitov's $\mathcal A(m,k)$, with Hermitian $A$ and at most $k$ nonreal $b_j$, is not preserved.
   - The same holds for the numerical normalization $t=b_1$ with nonreal $b_1$, and for the reciprocal transform with nonreal $c$.

   For example, the conditional route in [`../literature/README.md`](../literature/README.md) runs from a 4×4 Hermitian case to all $n$ via Shitov's Theorem 19. A proof on that route must establish **every** support direction θ directly, not only θ=0. This is a limitation on use, not an error in the reductions file.
2. **The rational-witness enumeration halts only on failure.** The file already says this. I repeat it because it is the easiest claim to over-cite: a bounded rational search that finds nothing is computational evidence only.

## 4. Spot checks

[`identity_spotchecks.py`](identity_spotchecks.py) tests the identities above in floating point: seed `20261125`, 200 trials, $n\in[2,6]$, tolerance $10^{-9}$. It covers:

- LUR invariance, spectrum swap, opposite shifts, and $c^n$ scaling.
- The reciprocal transform and its vertex scaling.
- The Cayley transform's unitarity, its inverse and the round trip.
- Kovačec's generalized-Givens determinant formula: $\det(D_a+GD_bG^*)=M\big((a_i+b_i)(a_j+b_j)+(a_i-a_j)(b_i-b_j)(d^2+s^2)\big)$.

The observed worst error was $4\times10^{-11}$, from the Cayley inverse near an eigenvalue $-1$. Every other worst error was at most $1.2\times10^{-14}$. These runs catch sign and index slips. They are not proofs.

```sh
python -I adversarial/identity_spotchecks.py   # numpy 2.4.6, Python 3.12.10
```

## 5. Not done in this session

- **The Zheng $n=4$ claim remains unaudited.** Its supplementary programs (S0–S2) are still not located; see the literature README. It should stay a claim pending audit. I suggest a separate task for it, because it is a large proof review rather than an assumption check.
- **The Corollary 3 subscript discrepancy is confirmed, not new.** The literature README reported it first. I rendered v2 p. 4 at 130 dpi and confirmed that the outside product is printed as $\prod_{\nu\ne i,j}(a_\nu+b_\nu)$ where direct expansion needs $b_{\pi(\nu)}$. I also confirmed the arithmetic of the literature's example: 343 against 441. A second slip on the same page: Corollary 3 is stated for a real Givens $G(i,j,c,s)$ but keeps the factor $d^2+s^2$, and $d$ is undefined there. With $d=0$ it should read $s^2$. Neither slip affects the segment Corollary 4 or MOC.
- **§5 of the reductions file and the minors / Theorem D artifacts were not reviewed.**
