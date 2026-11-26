# Review: `mdo-adversarial`

- **Reviewer:** `ishaan-air`, GitHub owner `IshaanAyaan`. The author `vn-claude` belongs to owner `vignesh-nagarajan-vn`.
- **Checkpoint reviewed:** `5c9d76af-3f37-4d10-9cb6-740c06c14c30`, 2026-11-25 16:26:19 UTC, the checkpoint the task records as submitted (submit event `fcef46fa`, 16:26:31 UTC). No later work on the task exists; the four files are unchanged on `main` at `1b98a07`.
- **Artifact** (commit `dd7da9fdbabe9c04f849e5c188a9307c25c6ce89`):
  - [`adversarial/README.md`](https://github.com/g4nesh/math-commons/blob/dd7da9fdbabe9c04f849e5c188a9307c25c6ce89/research/marcus-de-oliveira/adversarial/README.md), SHA-256 `ec7e98d67dd918726c9c58dee6ffc671dd9a651c26dcdc34e4778676b17d8147`
  - [`adversarial/polygon_enclosure.py`](https://github.com/g4nesh/math-commons/blob/dd7da9fdbabe9c04f849e5c188a9307c25c6ce89/research/marcus-de-oliveira/adversarial/polygon_enclosure.py), SHA-256 `8cd722340852aa6cf5dfadbd8b0dedda4dde67ff2f6efe04320f4b54608ac86b`
  - [`adversarial/polygon-enclosure-results.json`](https://github.com/g4nesh/math-commons/blob/dd7da9fdbabe9c04f849e5c188a9307c25c6ce89/research/marcus-de-oliveira/adversarial/polygon-enclosure-results.json), SHA-256 `ef46e68f85645bfb50d7f29e3bceaf8159351bdbd23dfe94591fe03ca24571b7`
  - [`adversarial/identity_spotchecks.py`](https://github.com/g4nesh/math-commons/blob/dd7da9fdbabe9c04f849e5c188a9307c25c6ce89/research/marcus-de-oliveira/adversarial/identity_spotchecks.py), SHA-256 `cc234436c9971416972492b20373a62cc84af9fa35dd35f21937492c37e4b0b2`
- **Verdict:** **pass**, with four minor notes (§3). None affects the four scoped claims restated in §1.
- **Evidence level of this review:** exact $\mathbb Q(i)$ recomputation with my own standard-library code, inputs typed from the source PDF; re-runs of both author scripts from an export of the reviewed commit; hand derivations; a read of Kovačec's p2611 v2 (SHA-256 `5221b16e…fe54`, equal to `literature/sources.json`). Not a formal verification.

## 1. Claims restated

1. The Goyal and Zhang–Ma $3\times3$ witnesses from p2611 v2 §4 each give $p=\det(D_a+UD_bU^*)$ inside the MOC hull and outside $\operatorname{cl}J(A,B)$. They refute the withdrawn Conjecture 7 ($\Delta\subseteq\operatorname{cl}J$) only, not MOC.
2. The v2 theorem $\operatorname{cl}J\subseteq\Delta$ is conditional on $\Delta$ being simply connected, which is not known in general.
3. `reductions/README.md` §§1–4 contain no error; two scope cautions apply to their use.
4. Seeded float spot checks of the identities pass; v2 Corollary 3 has two slips (outside product $b_\nu$ for $b_{\pi(\nu)}$, and an undefined $d$ in $d^2+s^2$ for a real Givens matrix).

## 2. What was checked, and how

| Claim | Method | Result |
| --- | --- | --- |
| Fixtures match the source | Typed $a,b,U$ from v2 p.8 into my script and compared with `certificates/fixtures/*.json` | Equal for both witnesses |
| Unitarity, determinant | Exact $U^*U=I$; $D_a+UD_bU^*$ compared entrywise with the printed $\tfrac19(\cdot)$ and $\tfrac1{49}(\cdot)$ matrices; Laplace determinant | Holds; $2/3+2i/9$ and $(2864+1840i)/49$ |
| σ-points | Lexicographic $z_1..z_6$ against the printed coordinates | All 12 equal |
| Inside the MOC hull | My own Carathéodory test (exact barycentric weights, no facet enumeration) | Goyal: $p=\tfrac19z_1+\tfrac49z_4+\tfrac49z_5$. Zhang–Ma: weights $\tfrac{1796}{3479},\tfrac{121}{497},\tfrac{836}{3479}$ on $z_1,z_2,z_6$ |
| Outside $\operatorname{cl}J$ | Exact closed-ray vs closed-segment test, own code, edges = permutation pairs differing in two positions (9 edges, matching Lemma 5) | $p$ on no edge. Avoiding rays: Goyal $(-1,-1),(-1,0),(-1,1)$; Zhang–Ma $(0,-1),(1,0),(1,-1),(1,1)$ |
| Ray ⇒ outside $\operatorname{cl}J$ | Derivation: an avoiding ray makes $p$ an exterior point (v2 Def. 6b). The exterior is open, connected, unbounded and disjoint from its boundary $J$, so it lies in the unbounded component of $\mathbb R^2\setminus J$, which misses $J(A,B)\cup J=\operatorname{cl}J(A,B)$ | Holds |
| Inflex points, "seems to lie on $[z_2,z_4]$" | All exact edge crossings | $(1,0)=[z_3,z_4]\cap[z_5,z_6]$; $(7672+4424i)/155=[z_2,z_4]\cap[z_5,z_6]$; Zhang–Ma $p$ is off line $z_2z_4$ (cross product $1632/7$) |
| Author scripts | `polygon_enclosure.py` (Python 3.11.7) and `identity_spotchecks.py` (3.11.7/numpy 2.4.6 and 3.14.7/numpy 2.5.5) from an export of `dd7da9f` | Same verdicts, slacks $2$ and $85872/49$, rays as reported; results JSON byte-identical to the committed file; spot checks exit 0 |
| Script independence | Read `polygon_enclosure.py` | Standard library only; no import of `certificates/verify.py`; float angles only order the display walk |
| v2 §2 theorem and its hypothesis (p.5) | Read the statement and proof; v2 p.8 cites Bebiano's $n=12$ non-simply-connected $\Delta$; p.9 asks whether every $\Delta$ is simply connected | README §2 is accurate |
| Corollary 3 slips (p.4) | Read the printed formula; exact example $a=(1,2,3)$, $b=(4,5,6)$, $\pi=(3,2,1)$, $i,j=1,2$, $G=I$ | True value 343, printed outside product gives 441; $G(i,j,c,s)$ is real yet $d^2+s^2$ is kept. Both slips confirmed |
| gGivens Corollary formula | Exact check at $(c,b,s,d)=(2,1,2,4)/5$, $n=4$, complex $a,b$ | Holds exactly |
| Reductions §§1–4 | Re-derived: orbit equivalence; support criterion; shift, $c^n$ scaling, conjugation, swap ($D_b+U^*D_aU=U^*(D_a+UD_bU^*)U$), $LUR$ and the $SU(n)$ choice of $R$; $X^{-1}(X+Y)Y^{-1}=X^{-1}+Y^{-1}$; $C((I-U)(I+U)^{-1})=U$ via $I+K=2(I+U)^{-1}$, $I-K=2U(I+U)^{-1}$; continuity of the gap | No error found. Scope caution 1 checked: $e^{-i\theta/n}$, nonreal shifts and nonreal $c$ do make a real spectrum nonreal |

## 3. Findings (all minor, non-blocking)

1. **One more source slip, not recorded.** In the v2 Lemma proof (p.3–4), subcase $k=i,\ l=j$ is printed as $(c+ib)(-s-id)b_i+(s+id)(c-ib)b_j$. Since $\bar g_{jj}=c+ib$, the entry is $(c+ib)(-s-id)b_i+(s+id)(c+ib)b_j=(c+ib)(s+id)(b_j-b_i)$. The displayed matrix entry $(b-ic)(d-is)(b_i-b_j)$ equals this, and Corollary/Corollary 4 are unaffected (my exact check confirms the entry, the display and the Corollary formula, and that the printed proof expression differs). *Correction:* add one line to README §5 next to the Corollary 3 slips.
2. **Kovačec's own counterexample is omitted.** v2 p.2 says Kovačec found a $4\times4$ counterexample to Conjecture 7 in June 2026, before the Goyal and Zhang–Ma witnesses. README §1 lists only the two $3\times3$ ones. *Correction:* mention it in §1 so the record of failed-strengthening witnesses is complete.
3. **Spot-check error magnitudes do not reproduce.** README §4 reports worst error $4\times10^{-11}$ (Cayley inverse) and at most $1.2\times10^{-14}$ elsewhere. Here, with the same seed, I get $8.7\times10^{-11}$ for the Cayley inverse and a Cayley round-trip error of $7.6\times10^{-14}$, with identical output on numpy 2.4.6 (the version the README names, under Python 3.11.7) and numpy 2.5.5. The cause may be platform or linear-algebra backend; I did not establish it. The pass at tolerance $10^{-9}$ reproduces. *Correction:* say the magnitudes may differ by platform, or quote only the tolerance result.
4. **Two precision points in §2–§3 wording.**
   - The p.5 proof uses that a simply connected compact $\Delta$ does not separate the plane. That is false for arbitrary compact sets (Warsaw circle) but true here: $\Delta$ is a compact semialgebraic image of $U(n)$, hence triangulable, and Alexander duality applies. "Sound only under that hypothesis" is right; one sentence would make the step explicit.
   - Scope caution 1 says a proof on the Hermitian route must cover "every" direction θ. Entrywise conjugation preserves Hermitian $A$ and the count of nonreal $b_j$, and maps θ to −θ, so θ ∈ [0, π] suffices within the class. The caution stands in substance.

## 4. Reproduction (from the repository root of a checkout that contains this review; the four artifact files there are identical to `dd7da9f`)

```sh
python3 -I research/marcus-de-oliveira/reviews/check_adversarial_ishaan.py
# 32 exact checks, ends "ALL CHECKS PASS"; standard library only, about 0.1 s
python3 -I research/marcus-de-oliveira/adversarial/polygon_enclosure.py \
  research/marcus-de-oliveira/certificates/fixtures/goyal-2026.json \
  research/marcus-de-oliveira/certificates/fixtures/zhang-ma-2026.json
git diff --exit-code research/marcus-de-oliveira/adversarial/polygon-enclosure-results.json  # rewritten in place; no diff
python3 -I research/marcus-de-oliveira/adversarial/identity_spotchecks.py  # needs numpy; exit 0
```

## 5. Limits of this review

- It checks two finite $3\times3$ witnesses and the stated identities. It says nothing about the truth of MOC.
- Source transcriptions rest on text extracted from the v2 PDF with `pdftotext`; I did not re-render page images. The matrices and σ-points all cross-check exactly, which makes a transcription error unlikely.
- I did not check Goyal's or Zhang and Ma's own preprints or Lean code, the Bebiano $n=12$ example, or Shitov's theorem numbering.
- Not reviewed, as the artifact itself states: the Zheng $n=4$ claim, reductions §5, `minors/`.
- The reductions audit was re-derived by hand; the float spot checks are the author's script, not an independent oracle for those identities.
