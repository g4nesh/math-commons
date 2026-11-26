Evidence level: EXACT for every identity, orbit membership, allocation and Gram coefficient relation listed below (own code, 211 checks, 0 FAIL); PROOF-READ for the PSD/closure steps that are pure linear algebra; NUMERIC (boundary attacks over SO(4)) as supporting evidence only; CITED where a step rests on part 2a (T24 face) or on part 1 atoms.

# SO(4) atoms, part 2b: allocation lemmas (G5) and the face Gram statements of historical T44, T45, T46, T47 (G1-G4)

Files: `so4_atoms_part2b.py` (own code; reads only supplement JSON data via `MOC4_SUPP`), `so4_atoms_part2b.json` (per-statement records and every exact check), `so4_atoms_part2b_output.txt` (captured run, 42 s).
Run: `MOC4_SUPP=<supplement root> MOC4_WRITE_JSON=1 python -I -B so4_atoms_part2b.py [G5 GRAM G1 G2 G3 G4 NUM]`.

## Deviations and limits

1. Scope is the universal SO(4) inequalities (atoms, allocations, Gram capacities) used by the listed documents. Not checked here, by design: the radial/polynomial parts of those documents (full-Hessian necessity G = 0, the positive-coefficient chart certificates, Bernstein certificates, KKT chart *coverage*, the "Gram-cone classification" and "exactly" claims of T46 core Thm 2 / six-ray). Those are UNCHECKED in this file.
2. Steps resting on the T24 face (part 2a) are CITED, not re-proved: T38/T27 eighth chamber (uses T24 coherent (5)), T38 gen.40 parallel/coherent capacities, T36 (7), T47 IIIb_minus_swapped (T24 negative-product coherent), T47 (3). Only the group transports into those roles were checked (EXACT).
3. Gram identities are checked in the coordinates the documents use, which the script rebuilds from supplement JSON *data* (cubic features, sign maps, the budget congruence T, and for T47 the native matrix V). Each identity is then verified exactly against the full linear system (rank 18, 26-dimensional solution space for every face), so the conclusion does not depend on where those polynomials came from. For T47 the documents' "spectral" coordinates had to be identified: they are the budget coordinates composed with V^-1 (g = V^-1 T^-1 F); the PCBE matrix satisfies all relations identically only there (other candidate congruences leave 15-18 relations violated).
4. Numeric attacks: 20 boundary samples per statement (coefficients exactly on the stated threshold), each minimised over SO(4) by 6000 random quaternion pairs plus 6 BFGS polishes, worst point re-evaluated in mpmath at 40 digits; seeds 1000+i, tolerance -1e-9 on Phi / max|coef|. This is sampling, not a cover of the region. Several thresholds are sufficient but not sharp (T44 Prop. 2 at b=j=c=1 survives p halved), so sensitivity controls use the sharp old-H atom lowered by 10% (found -7.3e-7) and the faces with their main minor removed (found -0.04 to -0.06).
5. The T45 literal old-H (row, col) tables were not reproduced by any of our readings of the action convention; the needed orbit membership with the stated roles holds (EXACT), so this is a notation point, not a gap.
6. In the coordinates used, minor f15 (T44) and f5 (T46) are not diagonal (|u_f|^2 is not the coordinate square); no check depends on it.
7. Process: one early read-only inspection was run outside the job limiter; every later run used it.
8. G4 items not separately checked: the old core c >= 1/e + k as a specialisation, PCBE_low_e_q_monotone's monotonicity argument, high_j_small_b conic certificates (not universal). Marked UNCHECKED below.

## G5. Allocation lemmas (T28, T38, T42, T36, T27, T48)

Statements and checks (all EXACT):
- Source identities: T28 (5), T38 (3), T42 (4), T36 (3), T27 (3) response identities hold on the whole contact nullspace; T36 contact nullspace has dimension 8, delta35 + delta40 = -g17, and decompositions (1), (2), (3) hold. T42 symmetry exchanging f22, f54 is realised (literal reading g = ((1,0,3,2),(0,2,1,3),0) matches).
- Transports: T24 core into T38, T27, T36; T28 face into T36; SOC into T36; the T36 stabiliser exchanging (a,b,q,p).
- H/chain allocation (6): with d = 1, beta1 + beta2 = 1, 0 <= delta_i <= beta_i, the cost b f_b + c f_c + a f_a + p f_p + q f_q + f_k is a sum of two old-H atoms and two chain atoms, each exactly at its threshold, whenever a >= 1 - d1 - d2 + b1^2/b + b2^2/c + d1^2/p + d2^2/q. Checked for the role tuples of T28 Lemma 2, T38 gen.70, T42, T36 (10); SOC partners in the proved orbit. KKT charts HH, HC, CC: energy equals the stated J, stationarity at beta1 = theta, both leaf slacks positive. Chart coverage: UNCHECKED (radial part).
- Single atoms: T28 (3) H(45,1;34), chain(1;0,45;34), four(6,0;5;34); T38 H(2,5;40), four(0,67;2;40), four(0,7;12;40); T27 four(67,41;39;40); T36 (6) SOC(45,1;35), H(5,45;35), four(7,0;5;35). All in proved orbits with stated roles.
- T27/T38 eighth chamber: the four-coefficient atom is exactly on its boundary; coherent scaling gives r^2/e - r^2/E_r; p(delta) increasing. Rests on T24 coherent (5): CITED (part 2a).
- T48 allocations I, II, III, IV: atoms H(2,5;40), H(29,5;40), chain(5;0,2;40), chain(5;7,29;40), chain(5;41,29;40), exact bookkeeping; Cor. 2 e-side identity holds.
Attack: H/chain (6) with T42 roles and the RHS minimised over (beta, delta): min 1.4e-18 (no violation).
Status: REDUCES TO PROVED ATOMS (allocations); PROVED HERE (identities, transports); eighth chamber REDUCES TO PROVED ATOMS + CITED part 2a.

## G1. Historical T44 (face q f2 + p f1 + j f53 + c f12 + e f8 + d f18 + b f7 + a f15)

- Prop. 1 (coherent chain family): for r != k, w > 0, with x = (v(2r-k)-1)/(r-k), y = k(1-rv)/(w(r-k)), the cost is nonnegative on SO(4) when p >= x^2+y^2, j >= r^2(1+v^2), c >= k^2(1+v^2), a, e >= v^2, b >= w^2, q >= 0, d = 1. Proof: explicit vectors J=(rv,r,0), C=(kv,k,0), E=A=(0,v,0), B=(0,0,w), P=(0,x,y), D=(1,-v,0); their Gram matrix satisfies every coefficient relation of the face (EXACT), so the cost is a sum of squares; slack in coefficients adds minors. PROVED HERE.
- Prop. 2: specialisation v=0, r=sqrt j, k=-z, w=sqrt b gives p(z) = (1+z^2/b)/(sqrt j+z)^2; derivative, both branches ((b+c)/(b(sqrt j+sqrt c)^2) if jc <= b^2, 1/(j+b) if jc >= b^2), j = 0 (1/b+1/c) and b = 0 (direct Gram, p = 1/j) all EXACT; the worked example reproduced. PROVED HERE.
- Scalar test (1): Prop. 1 at v = t/sqrt(1-t^2) etc. gives the displayed threshold. PROVED HERE.
- SOC(1,53;18) and Hsum (old-H (15,53;18) + (8,53;18), d split a:e): REDUCES TO PROVED ATOMS.
- c <= j note: transport (1,53,12,7,18) -> (53,1,8,7,18) realised in G (literal reading matches), so J* is the G-image of Prop. 2. PROVED HERE.
- native_contact_chart: no universal SO(4) inequality.
Attacks: Prop. 1 family min 2.9e-20; Prop. 2 threshold min 6.2e-19. No violation.

## G2. Historical T45 (contacts {3,4,5,9,10,12}; p f0 + r f1 + a f2 + b f12 + c f59 + q f67 + d f70 + s f72)

- (1): contact nullspace has dimension 8; the eight sources vanish on the contacts; their values at permutations (18,16,13,8,2,14,11,7) form the identity matrix; every nullspace cost equals the stated expansion; seven sources are minors, f70 is not. EXACT.
- (2)-(5): SOC(2,59;70), chain(2;0,12;70), chain(59;0,1;70), chain(59;67,72;70) in proved orbits with stated roles: REDUCES TO PROVED ATOMS.
- (6) c >= d + d^2/(r+s): H(1,59;70) and H(72,59;70) in the orbit with stated roles; splitting d in ratio r:s puts both atoms exactly at threshold and sums exactly to the target (EXACT). REDUCES TO PROVED ATOMS.
- Zero faces a=0, c=0, r=0: T33, T6, T14 contact sets are carried by some g onto the stated images. EXACT.
- Small PSD8 Grams of the note: not separately rebuilt (the atom reduction already proves (2)-(5)).
Attack: (6) split min 1.7e-18. No violation.

## G3. Historical T46 (face q f27 + p f1 + j f59 + c f5 + e f2 + d f28 + b f6 + a f4)

- Core coherent Thm 1: Gram G0 = [[q,beta,-t,0],[beta,b,-t,0],[-t,-t,j-t^2,0],[0,0,0,1]] in the basis (Q,B,J0,D) with P = lambda(Q+B) + mu J0, J = J0 + tD satisfies every coefficient relation (EXACT); |P|^2 equals threshold (1); G0 PD iff q>0, l = qb-beta^2 > 0, j > (1+(q+b-2beta)/l) t^2 (Schur, EXACT). Form (2), the optimal t*, u - 4/s = (q-b)^2/(s l), Cor. 3 (1/j + 2 - q - b), the beta = -1 forms (lambda, K, K - lambda, t = j/lambda, balanced case): EXACT. Finite limits by closedness of a fixed-diagonal slice: PROOF-READ. PROVED HERE.
- Six-ray Thm 1: G6 = [[q,beta,t,0],[beta,b,v,0],[t,v,j-t^2,0],[0,0,0,1]], P = lambda6 S + mu6 J0, C = kc(J0 - (v/b)B), D = C - Z: all relations, |P|^2 = p*, |C|^2 = c*, Schur condition (EXACT). PROVED HERE.
- Singular leaf (q != b): explicit PSD vectors, p* and c* formulas and the rational section p >= k + D0^2/(jN). PROVED HERE.
- four(4,6;2;28): REDUCES TO PROVED ATOM.
- Transports T24 -> T46 -> T47 and six-ray -> T47 [0,45,5,9,29,35,7,36]; six-ray Cor. 2 and Cor. 4 substitutions. EXACT.
- Core Thm 2 / six-ray "exactly" (cone classification): UNCHECKED (not a universal inequality).
Attacks: core min 1.1e-18; six-ray min 2.1e-18; control with p removed found -0.054.

## G4. Historical T47 (face q f2 + p f9 + j f1 + c f45 + e f5 + d f35 + b f0 + a f7)

- PCBE Thm 1 (d=1, p,b,j,e>0, 0<L<=min(p,b); c >= min_x max{K(1-x)^2-1, Bx^2}, K = 1/j+1/(L+1), B = max{1/e+1/(L+1), 1/L}): the PCBE matrix satisfies all coefficient relations identically in x, y, z, w, L (EXACT, spectral coordinates, see limit 3); rows P, B coincide and the Schur complement over {B,D} is the displayed 3x3 matrix (EXACT); scalar identity (3) (EXACT); simplification to K(1-x)^2-1 and Bx^2 (EXACT); x* solves the crossing (EXACT); (1) <=> (2) checked on 400 random triples (NUMERIC). The free JE entry completion (path graph J-C-E: PSD completion iff both 2x2 blocks PSD; explicit w given), admissibility in (3), clamping and closure: PROOF-READ, correct. PROVED IN THE DOCUMENT AND CHECKED.
- IIIb_plus_large_pc Lemma 1 (pc >= a^2, j >= d^2/(a+c)): Gram with A perp P,C, <P,C> = -a, J = d(C-A)/(a+c), D orthogonal satisfies all relations (spectral coordinates; the budget coordinates leave 16 violated) and is realised by explicit vectors in R^4 (EXACT). PROVED IN THE DOCUMENT AND CHECKED.
- merged-C Lemma 1: crossing, stationary and half-energy Grams for leaf a and leaf b satisfy all relations (6/6, EXACT); crossing value equals d^2(L+e)/(L(sqrt j+sqrt e)^2) (EXACT); PSD of each by the stated vector constructions (crossing: E parallel to -J; stationary needs je >= L^2; half-energy needs e >= j): PROOF-READ; proportional merge by degree-one homogeneity (EXACT). PROVED IN THE DOCUMENT AND CHECKED.
- Atoms used by I_II / high_j notes: H(9,1;35), SOC(45,1;35), chain(1;7,9;35), four(0,7;5;35) in proved orbits; transport T24 core (59,53,55,67,54,2) -> (45,5,7,0,35,1) realised. REDUCES TO PROVED ATOMS (+ CITED part 2a for the T24 role).
- zero_q_rational: reduces to six-ray Cor. 3 via the transport checked in G3. REDUCES TO PROVED ATOMS (G3).
- IIIb_minus_swapped and T47 (3): CITED (T24 negative-product coherent, part 2a).
- old core c >= 1/e + k specialisation, PCBE_low_e_q_monotone monotonicity, high_j_small_b conic: UNCHECKED.
Attacks: PCBE min 3.3e-18; large_pc min 5.2e-18 (control with j removed -0.056); merged-C min 8.1e-18. No violation.

## Summary of status

| item | status |
|---|---|
| G5 allocations (T28, T38, T42, T36, T48) | REDUCES TO PROVED ATOMS (EXACT bookkeeping) |
| G5 identities, transports | PROVED HERE |
| G5 eighth chamber | REDUCES TO PROVED ATOMS + CITED part 2a |
| G1 T44 Prop. 1, Prop. 2, scalar test, J* | PROVED HERE |
| G1 T44 SOC/Hsum | REDUCES TO PROVED ATOMS |
| G2 T45 (1)-(6), zero faces | REDUCES TO PROVED ATOMS / PROVED HERE |
| G3 T46 core, six-ray, singular leaf | PROVED HERE |
| G4 T47 PCBE, large_pc L1, merged-C L1 | PROVED IN THE DOCUMENT AND CHECKED |
| G4 T47 atoms, zero_q_rational | REDUCES TO PROVED ATOMS |
| G4 T47 IIIb_minus_swapped, (3) | CITED (part 2a) |
| radial / chart / classification claims | UNCHECKED (out of scope) |

No statement was found FALSE and no NOT PROVED item was found among the universal inequalities examined.
