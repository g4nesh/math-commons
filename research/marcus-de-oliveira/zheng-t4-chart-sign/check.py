#!/usr/bin/env python3
"""Independent exact replay of the T4 chart and scalar sign inputs.

No supplement programs, JSON certificates, or random trials are imported.
All equalities are identities of rational functions over QQ. Guards remain
active under python -O. This checks algebra; the inequalities are proved in
README.md, not inferred from a finite sample.
"""
from itertools import permutations
import json
from pathlib import Path
import sys
import sympy as s


def require(condition, label):
    if not condition:
        raise ValueError(label)


def equal(lhs, rhs, label):
    require(s.cancel(lhs-rhs) == 0, label)


def hermitian_matrix(poly, v, V):
    poly = s.cancel(poly)
    require(s.denom(poly) == 1, 'Hermitian polynomial has a denominator')
    p = s.Poly(poly, v, V)
    require(all(i <= 1 and j <= 1 for (i,j), _ in p.terms()), 'not bidegree (1,1)')
    require(p.coeff_monomial(v) == p.coeff_monomial(V), 'not Hermitian')
    return s.Matrix([[p.coeff_monomial(1), p.coeff_monomial(v)],
                     [p.coeff_monomial(V), p.coeff_monomial(v*V)]])


def replay():
    x,t,r,u,v,V,z,w = s.symbols('x t r u v V z w')
    h,y = 1-t,x-1
    H,L = h+u,y-r+u
    E,J = H*(1-t*x)-h*r*y,H*(1-t*x)+h*y
    denominator_z = r*t*v*v-r*v*x+t*t*v*v-t*u*v*v-t*v*x-t*v+u*v+x
    denominator_w = -r*t*v+r*x+t*u*v-u
    z_chart = r*v*(t*v-x)/denominator_z
    w_chart = (x-t*v)*(r-t*v+1)/denominator_w
    r_raw = z*(w-v)*(1-t*v)/(v*(z-w))
    u_raw = -(x-t*v)*(z*(1+v-w)-v)/(v*(z-w))
    raw = {r:r_raw, u:u_raw}
    dz_factor = -(v-w)*(x-t*v)*(t*v-1)/(z-w)
    dw_factor = -(z-v)*(x-t*v)*(t*v-1)/(v*(z-w))
    equal(denominator_z.subs(raw),dz_factor,'first chart denominator')
    equal(denominator_w.subs(raw),dw_factor,'second chart denominator')
    equal(z_chart.subs(raw),z,'chart recovers z')
    equal(w_chart.subs(raw),w,'chart recovers w')
    equal(r_raw.subs({z:z_chart,w:w_chart}),r,'chart recovers r')
    equal(u_raw.subs({z:z_chart,w:w_chart}),u,'chart recovers u')
    A=[z,0,1,x]; C=[w,v,None,t*v]
    products=[-s.prod(A[i]-C[p[i]] for i in range(4) if p[i]!=2)
              for p in permutations(range(4))]
    ratios=[s.cancel((p-products[0])/(v*(z-w))) for p in products]
    for j,expected in [(0,0),(1,1-x),(4,t-x),(5,x*(t-1)),(7,r_raw-x+1),(8,u_raw)]:
        equal(ratios[j],expected,f'raw contact {j}')
    ratios=[s.cancel(q.subs({z:z_chart,w:w_chart})) for q in ratios]
    gaps=[s.cancel((q-q.subs(v,V))/(v-V)) for q in ratios]
    require(len(gaps)==24,'all permutations')
    require(all(gaps[j]==0 for j in [0,1,4,5,7,8]),'all six contacts')
    a,b,c,d,e,k = gaps[15],gaps[19],gaps[10],s.cancel(gaps[2]-gaps[13]),gaps[14],gaps[16]
    identities=[(gaps[2],x/(v*V),'g2'),
                (gaps[6],r*t*y/((t*v-1)*(t*V-1)),'g6'),
                (c,-t*u*x/((t*v-x)*(t*V-x)),'c'),
                (gaps[9],y*c,'g9'),(a,y*c-d,'a'),
                (k,h*c-t*(b+d),'first source relation'),
                (E*e,t*y*H*b-J*d,'second source relation'),
                (J-E,h*y*(r+1),'J-E'),(J-x*h*H,-u*y,'J-xhH'),
                (E-H*(h-y)*(1+y),h*y*(H*y/h-r),'comparison E'),
                (L+u*(y-h)/h,H*y/h-r,'comparison L'),
                (E,(h-y)*(H+h*y)+h*y*L,'split E')]
    for left,right,label in identities: equal(left,right,label)
    Q=-t*(r+t-u)*v*v+(r*t+r*x+t*x+t-u)*v-x*(r+1)
    equal(Q,denominator_z*(z_chart-v)/v,'Q nonvanishing factor')
    normQ=Q*Q.subs(v,V)
    K=t*y*b-x*h*d
    MD=hermitian_matrix(d*normQ,v,V)
    MK=hermitian_matrix(K*normQ,v,V)
    equal(MD.det(),-r**3*t*t*u*x*y*L*E,'D determinant')
    equal(MK.det(),-r*r*t*t*u*u*x*x*h*y*L*E,'K determinant')
    q=s.Matrix([1,x/t])
    equal((q.T*MD*q)[0],r*u*x*y*(h*r-H*y),'D real evaluation')
    equal((q.T*MK*q)[0],r*u*u*x*x*y*(y-h),'K real evaluation')
    aa,bb,cc,dd=s.symbols('a b c d')
    capacity=aa*bb*(cc-dd)-dd*dd*(aa+bb+cc+dd)
    equal(capacity.subs(aa,y*cc-dd),cc*(y*bb*cc-x*dd*(bb+dd)),'capacity')
    arrow=s.Matrix([[aa+dd,0,-dd,0],[0,cc,dd,0],[-dd,dd,dd,dd],[0,0,dd,bb+dd]])
    equal(arrow.det(),dd*capacity,'arrow determinant')
    kk=s.symbols('k')
    slack=capacity.subs(aa,y*cc-dd)
    equal(slack-cc/h*((bb+dd)*(t*y*bb-x*h*dd)+y*bb*kk),
          y*bb*cc/h*(h*cc-t*(bb+dd)-kk),'quantitative slack identity')
    equal(dd-dd**2/(aa+dd)-dd**2/cc-dd**2/(bb+dd),
          dd*capacity/((aa+dd)*cc*(bb+dd)),'Schur complement')
    return {'status':'PASS','arithmetic':'exact rational functions over QQ',
            'sympy':s.__version__,'permutation_count':24,'contacts':[0,1,4,5,7,8],
            'denominator_z_factor':str(dz_factor),'denominator_w_factor':str(dw_factor),
            'D_matrix':[[str(s.factor(e)) for e in row] for row in MD.tolist()],
            'K_matrix':[[str(s.factor(e)) for e in row] for row in MK.tolist()],
            'scope':'Algebraic identities only; README contains sign proof. No global MDO verdict.'}


if __name__=='__main__':
    report=replay()
    if len(sys.argv)==3 and sys.argv[1]=='--output':
        Path(sys.argv[2]).write_text(json.dumps(report,indent=2)+'\n')
    elif len(sys.argv)!=1:
        raise SystemExit('usage: check.py [--output PATH]')
    print(json.dumps(report,indent=2))
