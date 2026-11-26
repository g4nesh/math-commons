# A complex second-entry rectangle with the fourth entry unrestricted

For the archived real four-cycle and fixed opposing spectrum, Marcus–de Oliveira holds for

$$A=\operatorname{diag}(x,\eta,y,z),\qquad x,y,z\in\mathbb C,$$

$$\boxed{\quad \frac14\le\Re\eta\le\frac34,\qquad |\Im\eta|\le\frac1{16384}.\quad}$$

This is a thin complex rectangle over the entire middle half of the real segment, including its boundary. The fourth eigenvalue remains unrestricted. Eighteen permutation masses are strictly positive; six stay zero. The exact evaluator returns all 24 masses.

The result extends the real segment of [PR61](https://github.com/g4nesh/math-commons/pull/61) into the complex plane. It complements the small square around zero in [PR59](https://github.com/g4nesh/math-commons/pull/59); neither region contains the other. It does not cover the whole second-entry plane or general $U,b$.

## Fixture and conventions

Use the same fixed data as [full-four-cycle](../full-four-cycle/README.md):

$$U=\frac1{625}\begin{pmatrix}-31&492&-240&300\\-300&225&500&0\\240&-180&225&500\\-492&-256&-180&225\end{pmatrix},\qquad b=(3-2i,2,1+2i,-i),\qquad B=U\operatorname{diag}(b)U^*.$$

The checker verifies unitarity. Thus $A$ and $B$ are normal. A permutation maps **columns to rows**, so $V_p(a)=\prod_j(a_{p(j)}+b_j)$. JSON uses lexicographic zero-based permutations. The six permutations with $p(4)=2$ are excluded. The [collision argument](../two-middle-positive-laws/README.md#fixture-convention-and-the-forced-face) forces them to have zero mass at $z=i$ when $\eta\ne i$; keeping them identically zero elsewhere is our certificate's choice.

Write $\eta=1/4+\tau/2+iv$, so $0\le\tau\le1$. Put $z=s+it$ and $D=1+s^2+t^2$.

## A uniformly positive real-interval law

[real-law.json](real-law.json) supplies 446 positive terms. Each is assigned to one of four cubic Bernstein controls, an allowed permutation, and a square $(a+bs+ct)^2$, with $a,b,c\in\{-1,0,1\}$. For each control and each of the eighteen allowed permutations, the coefficients on each of $1,s^2,t^2$ are at least one. Therefore their Bernstein combination gives

$$N_p^0(\tau,s,t)\ge D\quad(0\le\tau\le1).$$

The checker expands the determinant and pairing gaps using exact principal minors. It verifies every coefficient of the degree-four Bernstein gap identity. Thus this new strictly positive real law holds for arbitrary complex $x,y,z$ throughout $\Re\eta\in[1/4,3/4]$.

The earlier segment law was enough to prove its real interval, but some numerator Gram matrices were singular. The new strict numerator floor is what permits a uniform imaginary perturbation.

## Global rank control, proved by polynomial identities

Represent each allowed numerator by its six real quadratic coefficients in $s,t$. There are 108 coefficient variables. Expanding the four complex coefficients of $1,x,y,xy$ gives 80 real gap equations. Their matrix is

$$C(\tau,v)=C_0+\tau C_R+vC_I.$$

The row count includes zero coefficients; normalization is handled later by dividing by the sum of numerators. The checker reconstructs all three matrices from the original fixed $U,b$.

[rank.json](rank.json) gives **26 polynomial left-kernel relations**, each of degree at most three in $(\tau,v)$. Exact coefficient multiplication verifies $L(\tau,v)C(\tau,v)=0$. A recorded $26\times26$ minor of $L(1/2,0)$ is nonsingular. Hence $L$ has rank 26 over the rational-function field $\mathbb Q(\tau,v)$, so $C$ has rank at most 54 over that field. Every $55\times55$ minor of $C$ is consequently the zero polynomial. Therefore

$$\operatorname{rank}C(\tau,v)\le54$$

at **every real parameter pair**, even at points where the displayed relation matrix itself might lose rank. Independence is not assumed to persist from a single numerical sample.

## One minor stays invertible across the rectangle

Select the 54 rows and 54 columns in [bounds.json](bounds.json), and call their submatrix $B(\tau,v)$. Let

$$B_c=B(1/2,0),\quad X=B_c^{-1}B_R,\quad Y=B_c^{-1}B_I.$$

Exact elimination verifies $B_c$ is invertible. For $d=\tau-1/2$, we have $|d|\le1/2$ and

$$B(\tau,v)=B_c(I+dX+vY).$$

A finite geometric-series argument controls the whole real interval. The checker computes $I,X,\ldots,X^5$ and proves

$$T:=2^{-5}\|X^5\|_\infty<\frac{13}{50},$$

$$K:=\frac{\max_i\sum_{j=0}^4 2^{-j}\sum_k|(X^j)_{ik}|}{1-T}<20.$$

Indeed, factoring the inverse through $I-(-dX)^5$ gives $\|(I+dX)^{-1}\|_\infty\le K$ for every $|d|\le1/2$. This is a uniform analytic bound, not a grid test.

Let $c_0(\tau)$ be the 108 coefficients of the real law and let the subscript $R$ denote the selected rows. The four Bernstein controls of

$$V(\tau)=B_c^{-1}(C_Ic_0(\tau))_R$$

are checked exactly. They give $\|V(\tau)\|_\infty<314$ on the interval. The checker also proves $\|Y\|_\infty<28$.

At $r=1/16384$, the imaginary perturbation has contraction bound

$$r\cdot20\cdot28=\frac{35}{1024}<\frac12.$$

Thus the selected minor stays invertible for $|v|\le r$ throughout the full real interval.

## Correcting the law and retaining positivity

Adjust only the 54 selected numerator coefficients by

$$h(\tau,v)=-v(I+dX+vY)^{-1}V(\tau).$$

The selected 54 gap equations then vanish exactly. Since the full matrix has rank at most 54 and the selected minor is invertible, the selected rows span its row space: **all 80 gap equations vanish**.

The correction obeys

$$\|h\|_\infty\le\frac{|v|\,20\,314}{1-|v|\,20\,28}.$$

A real quadratic in $(s,t)$ with all six coefficients bounded by $H$ has a symmetric $3\times3$ Gram matrix with diagonal entries bounded by $H$ and off-diagonal entries bounded by $H/2$. Its operator norm is at most $2H$, so its absolute value is at most $2HD$. Consequently each numerator changes by at most

$$\frac{2r\cdot20\cdot314}{1-r\cdot20\cdot28}\,D
=\frac{785}{989}D<\frac45D.$$

Starting from $N_p^0\ge D$, every allowed numerator remains greater than $D/5$. The six excluded numerators remain zero. Set

$$H=\sum_pN_p>0,\qquad w_p=N_p/H.$$

These are nonnegative masses summing to one, and the full gap equations prove the asserted determinant identity. Individual quadratic coefficients after correction need not be positive; positivity of the **quadratic numerators** follows from the Gram bound.

The report also gives a uniform positive mass bound. If $E$ is an integer upper bound on the sum of square coefficients times squared linear-form norms at any real-law control, then $\sum_pN_p^0\le ED$. The correction contributes less than $72D/5$ in absolute total, so every allowed mass is at least $1/(5E+72)$. This is distinct from the numerator bound $D/5$.

## Replay and obtain weights

From the repository root, with Python 3.9+ standard library:

```sh
python3 -B research/marcus-de-oliveira/complex-second-region/check.py
python3 -B -m unittest discover -s research/marcus-de-oliveira/complex-second-region -p 'test_*.py'
```

The first command reconstructs the matrices, verifies the polynomial relations and exact independence, checks the finite-power and positivity bounds, evaluates original determinants, and compares [report.json](report.json). `--write-report` regenerates that report after the checks. Six tests cover certificate corruption, domain boundaries, exact input types, and off-grid original determinants with large complex fourth entries. `SHA256SUMS` pins the seven evidence files.

To evaluate the actual weights at $\eta=1/2+i/16384$, $z=100-100i$:

```python
import importlib.util
from fractions import Fraction
from pathlib import Path
path = Path('research/marcus-de-oliveira/complex-second-region/check.py')
spec = importlib.util.spec_from_file_location('rectangle', path)
rectangle = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rectangle)
print(rectangle.weights_at(Fraction(1,2), Fraction(1,16384), 100, -100))
```

The API requires exact integer or `Fraction` coordinates and rejects values outside the rectangle. The continuous proof applies to arbitrary real parameter values in it, not only rational API inputs. Discovery used numerical linear programming and FLINT, but neither is needed for replay.

## Attribution and limits

The fixed fixture and exact determinant arithmetic come from [four-cycle-spectral](../four-cycle-spectral/README.md). The real Bernstein construction follows [second-bernstein-segment](../second-bernstein-segment/README.md); this law adds a strict floor over a shorter interval. The coefficient and inverse methods build on [two-middle-positive-laws](../two-middle-positive-laws/README.md). The new contributions are the strict real law, polynomial global rank certificate, and one uniform rectangle with an actual evaluator.

On resume, authoritative main `b48a7ca` and Claude's [unitary factorization](../unitary-factorization/README.md) were checked. That result states a reduction to real orthogonal matrices and line/circle spectral corollaries. Our fixed opposing spectrum has three distinct collinear points and an off-line fourth, so it is neither collinear nor concyclic. Arbitrary $x,\eta,y,z$ need not lie on a line or circle either. Those corollaries do not establish this family and are not premises here. The archived $U$ is already real orthogonal.

Earlier floating rectangle attempts failed or timed out; they are not mathematical obstructions. The present evidence is a written exact proof with standard-library replay and same-owner tests, pending independent-owner review. It is not a Lean proof. The rectangle is deliberately narrow. The full complex second-entry plane, general $U,b$, and the general conjecture remain unresolved. A next extension is a larger analytic region using another well-conditioned minor or a less conservative Gram bound, with exact overlap and positivity checks.
