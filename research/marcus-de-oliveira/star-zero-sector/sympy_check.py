"""Optional independent expansion replay (requires SymPy; no numerical search)."""
import json
import sympy as S

V, W, Z, t, s, K0, R0 = S.symbols('V W Z t s K0 R0', real=True)
Q, K1, m = W+Z, V+W, 1-s
A, k, R, d = V*s, s*K1+m*K0, m*R0, Z*m

original = (A*A+A*(k+2*A+R)*t+(k*(2*A+R)+A*(d-W))*t*t
            +(k*d+Z*(R-m*K0))*t**3)
F0 = t*t*R0*(K0+Z*t)
F1 = V*V+V*(3*V+W)*t+V*(2*V+W)*t*t
Fc = (V*(K0+R0)*t+(2*V*K0+K1*R0+V*(Z-W))*t*t
      +Z*(R0+K1-K0)*t**3)
assert S.expand(original-((1-s)**2*F0+s*(1-s)*Fc+s*s*F1)) == 0

L = V+K1*t+Z*t*t
m11 = V*V*W+V*(V*W+Z*Z)*t+Z*(V*W+Q*K1)*t*t
m22 = V*(Z+Q*Q)+(Z*K1+Q*V*(W+3*Z))*t+V*Z*(Q+Z)*t*t
# Independent form assembly: extract the coefficients of K0, R0 and 1
# from Fc/t. K0=Q|G|² and Q R0=VW|B|²+Z|G|²+2vwz Re(conj(B)G).
constant = S.expand(Fc/t).subs({K0: 0, R0: 0})
k_coefficient = S.diff(Fc/t, K0)
r_coefficient = S.diff(Fc/t, R0)
raw11 = V*W*r_coefficient+Q*constant
raw22 = Z*r_coefficient+Q*Q*k_coefficient+Q*constant
assert S.expand(raw11-m11) == 0
assert S.expand((raw22-m22).subs(V, 1-W-Z)) == 0
assert S.expand(r_coefficient-L) == 0

certificate = (V*W*(V-Z*t)**2
    +(2*V**3*W+V**3*Z+V**2*W*Z+2*V**2*Z**2)*t
    +(V**3*W+2*V**3*Z+6*V**2*Z**2+2*V*W*Z**2)*t**2
    +(V**3*Z+2*V**2*W*Z+6*V**2*Z**2+V*W**2*Z+5*V*W*Z**2+W**2*Z**2)*t**3
    +(2*V**2*Z**2+V*W*Z**2)*t**4)
determinant = m11*m22-V*W*Z*L*L
assert S.expand((determinant-Q*Q*certificate).subs(V, 1-W-Z)) == 0
remainder = S.Poly(certificate-V*W*(V-Z*t)**2, V, W, Z, t)
assert all(coefficient >= 0 for coefficient in remainder.coeffs())
print(json.dumps({'sympy': S.__version__, 'bernstein_identity': True,
                  'independent_form_assembly': True, 'determinant_square_identity': True,
                  'remaining_coefficients_nonnegative': True}, indent=2))
