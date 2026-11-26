# When a single permutation law stops working for Fourier bases

For each ordinary normalized discrete Fourier matrix Fₙ with **n ≤ 8**, one permutation law reproduces every squared minor. Thus both eigenvalue lists in the Marcus–de Oliveira determinant can vary arbitrarily at each of these bases. **F₉ already fails at order two.** In fact the same pair constraint fails for every unitary within entrywise distance 1/4096 of F₉. The original determinant question still permits different weights for different spectra in that neighborhood.

The cases n ≤ 3 follow from the singleton law and complementary minors. Orders [five](../fourier-minor-laws/README.md) and [seven](../fourier-seven/README.md) have separate positive laws that persist on full neighborhoods. Below are exact constructions for four, six and eight, followed by the obstruction at nine. The construction concerns this particular list of Fourier bases. Obstructions for other unitaries already occur in smaller dimensions. Background on the stronger common-minor question is given by [Drury](https://www.sciencedirect.com/science/article/pii/002437959490362X). The literature priority of these Fourier-specific results is open.

## Sampling the positive laws

Let p send column indices to row indices in ℤ/nℤ. In each row of the table, use a seed with its indicated probability. Given p, draw h uniformly from the units modulo n, e uniformly from {1,−1}, and b,c independently uniformly modulo n. Output L p R⁻¹, where L(i)=hi+b and R(j)=e h⁻¹ j+c. This averages over relabelings that preserve Fourier squared minors. Tuple entries give p(0),…,p(n−1).

| n | Seed | Probability |
| --- | --- | --- |
| 4 | (3,2,0,1) | 1 |
| 6 | (5,4,1,0,3,2) | 1/4 |
| 6 | (5,4,0,2,1,3) | 1/8 |
| 6 | (5,3,1,4,0,2) | 5/8 |
| 8 | (7,6,5,1,0,3,2,4) | (3−2√2)/4 |
| 8 | (7,6,4,0,2,5,1,3) | √2−1 |
| 8 | (7,6,4,0,2,1,5,3) | (3−2√2)/4 |
| 8 | (7,5,2,0,4,6,1,3) | 1/2 |

Each list has positive masses summing to one. To verify the marginal for any I,J, count the relabelings for which L p R⁻¹ maps J to I and take the indicated weighted average. Independently expand the Fourier subdeterminant into signed powers of exp(2πi/n), multiply by its conjugate and divide by n to the power |I|. The [exact finite proof](check.py) performs both expansions for all 70, 924 and 12,870 pairs respectively, at every order including zero and n. All numerical coefficients belong to ℚ(√2) and the program checks the two rational coefficients separately. Numerical optimization selected candidate seeds. The verification uses no optimization and no rounding. The determinant conclusion follows by the [coefficient expansion](../fourier-minor-laws/README.md#why-the-minor-identities-give-the-determinant).

## A two-point statistic at order nine

For a permutation p of ℤ/9ℤ and an unordered pair {j,k} define d=(j−k)(p(j)−p(k)) modulo 9. Let T₂ count pairs with d=±2 and T₄ count pairs with d=±4. Let T₀ count pairs with d=0. The following linear inequality holds for every permutation and hence its expectation under every probability law:

    T₄ − T₂ − (3/2)T₀ ≤ 6.

An exact finite check is short. Subtract p(0) from every output, which preserves all d's and fixes p(0)=0. For each of the 8! permutations of the remaining eight elements, count d in the 36 pairs. These are the attained maxima:

| T₀ | 0 | 2 | 3 | 5 | 9 |
| --- | --- | --- | --- | --- | --- |
| max(T₄−T₂) | 6 | 9 | 10 | 8 | 9 |

For Fourier entries, the proposed probability of mapping {j,k} to {r,s} is (2−2 cos(2π(j−k)(r−s)/9))/81. When the product is zero this vanishes. For each of ±2 and ±4 there are 243 ordered choices of a source pair and a destination pair with product in that sign class. There are three unit difference classes up to sign, nine pairs of each class, and exactly one matching destination class for each source class. The 81 zero-product choices use the nine pairs of difference ±3 on each side. Hence, for α=2 cos(2π/9),

    E_F(T₀) = 0,
    E_F(T₄−T₂) = 3(2α²+α−4).

Here 2cos(4π/9)=α²−2 and 2cos(8π/9)=−α²−α+2, using α³−3α+1=0. The excess over six equals **3(2α−3)(α+2)**. The cubic is strictly increasing on (1,2), contains α there as its root, and is negative at 61/40. The excess therefore exceeds 1/2. This contradicts the permutation inequality using just the proposed two-element probabilities.

## Moving away from F₉

Let δ=max |U[r,j]−F₉[r,j]| ≤ 1/4096 for a unitary U. All entries of both matrices have magnitude below 1/2. A 2 × 2 determinant is a difference of two products and its change is at most 2δ. Both determinant magnitudes are at most 1/2, so the squared modulus changes by at most 2δ. The entire statistic therefore changes by at most

    2δ · (243 + 243 + (3/2)·81) = 1215δ < 1/2.

It still exceeds six throughout this ball. This proves an open set where even the order-two permutation marginals are incompatible. To advance the determinant conjecture there, a proof must allow spectrum-dependent weights or use a different sufficient condition.

Run `python3 -B research/marcus-de-oliveira/fourier-composite-laws/check.py` from the repository root to replay the positive laws, the exact finite inequality and the rational comparisons in the neighborhood bound.
