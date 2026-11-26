# The full U11=0 boundary of the four-dimensional inequality

**Claim for independent review.** The inequality ★ holds for every complex
unitary U of order four with U11=0. This does not prove ★ for arbitrary U,
or the Marcus–de Oliveira conjecture. The proof below is algebraic and
includes all phases, row completions and zero cases. No novelty claim or
Lean formalization is made.

The task began with two ordering sectors, hence the directory name. A
second two-by-two matrix argument removes the remaining ordering. The
only input is the exact scalar reduction in
[the two-row report](../star-two-row/README.md), pinned at
`36742932d10ee25012b11714678468917031f3f8`. This report supplies its own
polynomial identity replay; the matrix/compression bridge is written there.

## 1. Target and reduction

For one-based indices, let B_ij=|U_ij|² and D[S,T]=|det U[S,T]|². Put

    A=D[12,12],  C=D[13,14],
    g0=B11+B33+D[14,34]−D[12,14].

The target ★ is g0+A+2√(AC)≥0, equivalently
f_t/t=g0+(1+1/t)A+tC≥0 for every t>0. This includes A=0 or C=0 by taking
the infimum over t; division by either minor is never used.

Column phases make the first row x=(0,v,w,z) real and nonnegative. Set

    V=v², W=w², Z=z², Q=W+Z, K1=V+W;   V+W+Z=1.

If Q=0 then z=0. The original formula is
g0=|b3|²≥0, so ★ holds. Hence assume Q>0, put q=√Q, and write row two as

    a=α e1+β(0,0,z,−w)/q+γ(0,q,−vw/q,−vz/q),
    |α|²+|β|²+|γ|²=1.

Let s=|α|² and r=|γ|². The two-row report proves that positivity for all
remaining row completions is equivalent to a nonnegative trace and the
following cubic F(t):

    A=Vs, k=(1−Z)s+Qr, R=|vwβ+zγ|²/Q, d=Z(1−s),

    T(t)=2A+(k+2A+R)t+dt²,
    F(t)=A²+A(k+2A+R)t
          +[k(2A+R)+A(d−W)]t²
          +[kd+Z(R−Qr)]t³.

The trace T is nonnegative immediately. It remains to prove F(t)≥0.

## 2. Decompose in the occupation s

If s<1, write β=√(1−s)B and γ=√(1−s)G, where B,G are complex and
|B|²+|G|²=1. Define

    K0=Q|G|²,  R0=|vwB+zG|²/Q.

Then k=sK1+(1−s)K0 and R=(1−s)R0. Direct expansion gives the exact
Bernstein decomposition

    F(t)=(1−s)² F0(t)+s(1−s) Fc(t)+s² F1(t),

    F0(t)=t²R0(K0+Zt),
    F1(t)=V²+V(3V+W)t+V(2V+W)t²,

    Fc(t)=V(K0+R0)t
           +[2VK0+K1R0+V(Z−W)]t²
           +Z(R0+K1−K0)t³.

The two endpoint polynomials F0,F1 are nonnegative. When s=1, β=γ=0
and F=F1 directly, so no normalization of B,G is required at that endpoint.

## 3. The middle polynomial is a positive quadratic form

For t>0, define L=V+K1t+Zt² and the real symmetric matrix M(t) by

    M11=V²W+V(VW+Z²)t+Z(VW+QK1)t²,

    M22=V(Z+Q²)+[ZK1+QV(W+3Z)]t+VZ(Q+Z)t²,

    M12=M21=vwz L.

Substitute K0,R0 into Fc and use |B|²+|G|²=1 and V+W+Z=1. The result is

    Q Fc(t)/t = (B,G)* M(t) (B,G).

Both diagonal entries of M have only nonnegative coefficients. Its
determinant has the following exact square decomposition:

    det M = Q² {
      VW(V−Zt)²
      +[2V³W+V³Z+V²WZ+2V²Z²]t
      +[V³W+2V³Z+6V²Z²+2VWZ²]t²
      +[V³Z+2V²WZ+6V²Z²+VW²Z+5VWZ²+W²Z²]t³
      +[2V²Z²+VWZ²]t⁴
    }.

Every displayed term is nonnegative for V,W,Z,t≥0. Thus M is positive
semidefinite, including zero diagonal entries and zero determinant. Its
quadratic form is nonnegative for complex B,G as well as real ones.
Consequently Fc≥0, then F≥0, then the original compressed matrix is
positive semidefinite. This proves f_t≥0 for every t>0 and every remaining
row completion. Taking the infimum proves ★ on the entire U11=0 stratum.

No restriction on the ordering of V,W,Z, on the phases of α,β,γ, or on
whether a minor vanishes is used. Q=0 and s=1 were handled before division;
V=0, W=0, Z=0, s=0 and repeated eigenvalues of either two-by-two matrix
are all allowed by positive semidefiniteness.

## 4. Exact verification and provenance

From the repository root, using Python 3 and no external dependencies:

```sh
python3 -B research/marcus-de-oliveira/star-zero-sector/verify.py
```

The checker uses sparse polynomials with exact integer coefficients to
verify the Bernstein identity, the independent quadratic-form assembly,
and the determinant square identity modulo V+W+Z=1. These are symbolic
polynomial identities, not numerical samples. It also replays the earlier
complete Gaussian-rational U11=0 negative-coefficient matrix against the
new decomposition, and checks representative exact boundary matrices.
`--write-report` rewrites only the adjacent report file.

As a separately implemented algebra check, SymPy 1.14.0 verified all three
identities by expansion and substitution in an isolated local environment.
The public `sympy_check.py` reproduces that optional check if SymPy is
already installed; it is not needed by the standard-library checker.

The earlier negative t² coefficient is genuine: coefficientwise positivity
of F is false. This proof uses the Bernstein decomposition and the square
inside det M, which preserve the cancellation. No numerical search was
extended and no outside violation or universal ★ claim is inferred.

The derivation is by ganesh-reductions; the coordinator supplied the
two-row compression and requested the boundary follow-up. Same-owner
symbolic or algebra checks are not independent mathematical review.
