# Source check of reduction.md against CLMP arXiv 2505.24326v3 (raw LaTeX)

Download (2026-11-27 19:18 UTC), untrusted data, only read with grep/Read, nothing executed:

    curl -sL -o src.tgz https://arxiv.org/e-print/2505.24326v3 ; tar -xzf src.tgz -C x
    sha256(src.tgz) = e70e507f60e2b222c83e9680c541fb8a0f9963071d29a04e2e3a532cf8648310   (32425 bytes, gzip, mtime 2025-12-09)

Local dir: `C:\Users\senth\AppData\Local\Temp\claude\0833dcd534d9398f\1ea7be13-f9cc-42c5-bce6-303c8d4f39eb\scratchpad\arxiv-src\src1\x\` (main.tex, Sections/*.tex). Files below are relative to it. LaTeX labels are given because the compiled numbers (1.4, 3.6, 5.12) cannot be confirmed without compiling; only the label names are certain. The HTML fallback was not needed.

## a. Theorem 1.4 (label thm:inductionstep), Sections/IntroGP.tex:50-53
> Let $N=pN'$ be a square-free number with $p \in \NN$ prime and $N' \in \NN$ such that all principal minors of $\mc{F}_{N'}$ are non-zero in characteristic $p$. Then, any $N'$-principal minor (and thus any principal minor) of $\mc{F}_N$ is non-zero.

Hypotheses: N square-free (so p does not divide N'), p prime, all principal minors of F_{N'} nonzero in char p. No condition p vs N' (no p < N', no parity or size condition), no extra condition on N. Draft sec. 1/open point 3: **CONFIRMS**.

## b. N'-principal (Definition defn:minors), Sections/Fourier_minors.tex:13-26
> $K_{d,i}=\set{k\in K: k\equiv i\bmod d}$ ... The matrix $B$ is a \emph{principal} submatrix of $A$ if $K = L$; a \emph{$d$-principal} submatrix of $A$ if $\abs{K_{d,i}}=\abs{L_{d,i}}$ for all $0\leq i\leq d-1$
> It is clear that every principal submatrix is $d$-principal for any divisor $d$ of $N$

Since K=L gives equal counts in every class, every principal minor of F_N (any S, any size, incl. empty by convention) is N'-principal, and Thm 1.4 says it is nonzero. Also IntroGP.tex:61: "principal minors are a special case of $N'$-principal minors". Draft: **CONFIRMS** (draft says Def 3.2; the label is defn:minors, number not checked).

## c. "Nonzero in characteristic p" (Section 3.2, subsec:basic-facts-cyclotomic-fields; q-Chebotarev Definition), Sections/Number_Theory.tex
- 52-55: for q not dividing n, f=ord_n(q), "The Fourier matrix $\mc{F}_n$ in characteristic $q$ can be obtained by taking the entries of $\mc{F}_n$ modulo $\mf{Q}$, where $\mf{Q}=\mf{Q}_j$ with any $1\leq j\leq r$".
- 59: "$\ol{\mc{F}}_n := (\ze^{ij})_{0\leq i,j\leq n-1}$" for a fixed n-th root of unity zeta in F_{q^f}. So it is the **unnormalized** matrix (zeta^{ij}); no 1/sqrt(n).
- 70: $\ol{D}_{K,L}=\det (\ze^{ij})_{i \in K, \, j \in L} \in \FF_{q^f}$.
- 132: "if $\ol{D}_{K, L}=0$, then $D_{K,L}\in\mf{Q}$, and the Galois group ... acts transitively on the prime ideals dividing $q$, so that $q\mid \Norm(D_{K, L})$. Conversely, if $q\mid\Norm(D_{K, L})$, then at least one of the Galois conjugates of $D_{K, L}$ belongs to $\mf{Q}$."
- 184-193 (q-Chebotarev property, four conditions "equivalent" as the paper states): (1) "All minors of $\mc{F}_M$ are non-zero in characteristic $q$." (2) "All minors of $\ol{\mc{F}}_M$ (with an $M$-th root of unity $\ze$ chosen from $\FF_{q^f}$) are non-zero." (3) "...$q\nmid \Norm(D_{K,L})$." (4) "$\ol{D}_{K, L}\neq0$".

Verdict: **CONFIRMS with one nuance.** Theorem 1.4 never fixes a prime, and the proof (see d) uses the norm condition, i.e. nonzero modulo every prime above p. The paper's (3)<=>(4) equivalence holds because, on the full set of minors, the Galois conjugate D_{K,L} -> D_{K',L'} is again a minor, so "for one zeta, all minors" equals "for all primes". That argument at line 132 produces conjugates D_{kK,lL}, which are NOT principal in general. So for principal-only minors, the equivalence of "one prime" and "all primes" is not given by the paper. Checking both zeta and zeta^-1 (one representative per Frobenius orbit) is the correct, safe reading, as the draft does. Note the paper's condition (3)/(4) are stated for all minors (K,L arbitrary), whereas Thm 1.4 says only "principal minors ... non-zero in characteristic p". The proof (d) supports the principal-minor reading, using the norm.

## d. Proof of Thm 1.4, Sections/squarefree.tex:11-173; label eq:Nfinalevaluation at 162-167
Which finite-field minors are used (squarefree.tex:155-157):
> $\prod_{k=0}^{N'-1}(\det Z_k)^{M_{\sig(k)}-M_{\sig(k+1)}}\det V_{\sig(k),\sig(k)}$, where $M_{\sig(N')}=0$, and $Z_k$ is the square matrix formed by the first $k+1$ rows and columns of $Z=\bra{\ze^{\sig(i)\sig(j)}}_{0\leq i,j\leq N'-1}$. All $Z_k$ are principal minors of $\mc{F}_{N'}$ with permuted rows and columns

with $\ze=\om^p$ (line 50) a primitive N'-th root, and sigma any permutation sorting the class sizes M_a decreasing (129-130). So Z_k = F_{N'}[S_k,S_k] up to simultaneous reordering, S_k = {sigma(0),...,sigma(k)}: plain principal minors of the unnormalized (zeta^{ij}) matrix, not twisted. (Principal because the minor is N'-principal: rows and columns have the same class counts.)

Where the hypothesis enters (172): "Taking the norm $\Norm=\Norm^{\QQ(\ze)}_\QQ$ on both sides of \eqref{eq:Nfinalevaluation}, we get an integer divisible by $p$ on the left-hand side, whereas on the right-hand side we obtain an integer that is not divisible by $p$; indeed, $p\nmid\Norm{\det(Z_k)}$, $0\leq k\leq N'-1$ by hypothesis, $\Norm(\ze)=\pm1$ ..., $\Norm(\ze^a-\ze^b)\mid N'$ ..., each $m_{\sig(k),j}-m_{\sig(k),i}$ is not divisible by $p$, and the product of factorials is also not divisible by $p$."

So the proof needs p not dividing Norm_{Q(zeta_{N'})/Q}(det Z_k) for each Z_k that occurs. The norm is over all Galois conjugates, so the needed condition is nonzero modulo EVERY prime above p. Left side: Phi_N(zeta)=pu (Number_Theory.tex:225-229, eq:valuecyclotomic). Both zeta and zeta^-1 reductions (the two Frobenius orbits of primes above 11 in Q(zeta_14)) are needed; the proof uses a norm, so no choice of single prime suffices in general. Draft sec. 1 proof skeleton and open point 1: **CONFIRMS**. Also **CONFIRMS** the draft's worry that Z_k are plain (not Galois-twisted) principal minors. Equation label: eq:Nfinalevaluation is the 11th numbered display in this section by count (could be (5.11) or (5.12) depending on numbering; cannot confirm without compiling), same content.

Slack the draft may use: the proof only ever invokes the nested chain S_0 subset S_1 subset ... (prefixes of one ordering), not every principal minor, but the theorem's stated hypothesis is all principal minors of F_{N'}. Checking all 2^14 is sufficient and matches the theorem statement. Fine.

## e. Does the paper tabulate or state which small N' satisfy the hypothesis?
Searched all .tex (excluding comments): no mention of N'=14 or p=11 with N'=14, no active table. A table input is commented out (IntroGP.tex:93-99, `%\input{Sections/matrix4}`, `%\caption{... table indicates by 1 if q is p-primitive ...}`; matrix4 is not in the archive). Related active statements only: Prop N-6p (squarefree.tex:228-247) checks F_6 in char p>3; the stated resolved small cases: Theorem smallmultp (IntroGP.tex:43) "$N$ ... of the form $2p,3p,5p,6p$ or $7p$" and IntroGP.tex:46 "the smallest square free $N$ for which Conjecture ... is not resolved positively is, for two factors, $11 \times 13=143$, for three factors, $2 \times 5 \times 7=70$". So the paper (v3) does not claim N=154 or the (11,14) condition, and it lists N=70 as unresolved, which is the three-factor analogue. **UNRESOLVED** by the paper: the char-11 condition for F_14 must be computed. (The draft's remark that Theorem 1.3 covers F_14 over C is consistent with 2p, p=7; unaffected.)

## Effect on the computation
No change: check det(zeta^{ij})[S,S] != 0 in F_1331 for all 2^14 subsets S of Z_14 and for both zeta and zeta^-1 (the two primes above 11). Reminder that line 132 shows only the combination of all-prime with all-minor symmetry; for principal minors the draft correctly does not assume a single prime suffices. If any S fails for either root, this route fails (not the conjecture). Sources of risk: none found in the proof; the Section 5 proof is internally consistent on the points checked (Z_k principal, norm over Q(zeta_{N'}), m differences < N congruent mod N' so not divisible by p).
