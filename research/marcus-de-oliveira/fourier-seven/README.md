# Four permutation orbits certify the order-seven Fourier basis

Let F be the normalized 7 × 7 Fourier matrix with entries ζ^(ij)/√7 and ζ = exp(2πi/7). Both complex spectra can vary freely when the relative unitary U satisfies max |U[i,j]−F[i,j]| ≤ **1/32768**. The construction gives one probability law on all 5040 permutations depending only on U, with each mass greater than 1/161280. For every pair I,J of equal size the event p(J)=I has probability |det U[I,J]|². The determinant identity for every pair of spectra then follows by the principal-minor expansion and Cauchy–Binet as in the [order-five proof](../fourier-minor-laws/README.md).

## The law at the center

Use the labels 0 through 6 and do index arithmetic modulo seven. Put α = 2 cos(2π/7). Give the uniform law on all permutations mass 1/16. Use the remaining mass on the following four seeds. Each tuple lists the images of the seven column indices in order.

| Seed | Total mass of its randomized images |
| --- | --- |
| (6,5,4,3,2,1,0) | 1/40 |
| (6,5,4,1,0,3,2) | 137/48 − α − α² |
| (6,5,2,0,4,1,3) | −2309/240 + 3α + 4α² |
| (6,4,1,5,3,0,2) | 1843/240 − 2α − 3α² |

For each seed p, sample L(i)=hi+b and R(j)=e h⁻¹j+c, with h uniform on the six nonzero residues, e uniform on {1,−1}, and b,c independent uniform residues. Use the permutation L p R⁻¹. This averages over 588 transformations, counting repetitions. Translating row and column indices rephases the Fourier matrix, and opposite dilations preserve its entries up to conjugation. Thus each transformed pair of index sets has the same Fourier squared minor. The four masses add to 15/16. For positivity, α is the unique positive zero of x³+x²−2x−1 and lies strictly between 1.2469 and 1.2470. Bounding each polynomial term on that rational interval gives a positive lower bound for every displayed mass.

## The finite coefficient proof

The table gives enough integer counts to check the law by hand. A string such as 013 denotes the subset {0,1,3}. Entries in the middle column count, out of 588, the transformed versions of each of the four seeds with p(J)=I. The last column gives 7^k times the required Fourier squared minor at order k. For each row, multiply the counts by the four masses, divide by 588, and add 1/(16·binom(7,k)) from the uniform component.

| I, J | Counts for the four seeds | Scaled squared minor |
| --- | --- | --- |
| ∅, ∅ | 588,588,588,588 | 1 |
| 0, 0 | 84,84,84,84 | 1 |
| 01, 01 | 28,48,12,0 | 2−α |
| 01, 02 | 28,0,24,36 | 4−α² |
| 01, 03 | 28,36,48,48 | 1+α+α² |
| 012, 012 | 28,12,0,0 | 21−7α−7α² |
| 012, 013 | 0,12,12,12 | 7 |
| 012, 014 | 28,48,24,12 | 7α² |
| 012, 024 | 28,0,36,48 | 14+7α |
| 013, 012 | 0,12,12,12 | 7 |
| 013, 013 | 42,24,24,24 | 14 |

These are the representatives for the stated action on equally sized pairs of subsets through order three. Jacobi's identity gives the complementary orders. An independent verification dispenses with both the orbit classification and complementation. Run `python3 -B research/marcus-de-oliveira/fourier-seven/check.py`. The program computes the determinant by its signed permutation expansion in ζ, multiplies by the conjugate, and compares it with the four orbit probabilities for each of the **3432** pairs at orders zero through seven. In the real cubic field, ζ+ζ⁻¹=α, ζ²+ζ⁻²=α²−2, and ζ³+ζ⁻³=1−α−α². The checks compare all three rational coefficients exactly. Numerical linear programming suggested the four seeds. The proof and the replay use only their displayed values and exact arithmetic.

## A full neighborhood

Apply the signed determinant mass correction from the [order-five proof](../fourier-minor-laws/README.md#moving-the-basis). If δ is the distance from F in maximum entry norm, all entries of U and F have magnitude below 1/2 on the stated ball. Telescoping a seven-entry product bounds its change by (7/64)δ. Replacing columns one at a time bounds the change in the determinant by 7√7 δ. Multiplying that last bound by the unchanged Fourier monomial magnitude 7^(−7/2) contributes δ/49. Thus each signed mass changes by at most (7/64+1/49)δ. The reserve at F is at least 1/(16·5040)=1/80640, and at δ=1/32768 the possible change is less than 1/161280. The corrected law has the required positive masses throughout the closed ball.

This result covers the stated dimension-seven relative bases and their row and column rephasings and permutations. The general Marcus–de Oliveira problem allows arbitrary relative bases. Independent mathematical review and the literature priority of this construction remain open.
