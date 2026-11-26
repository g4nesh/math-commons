"""Exact spectral triangle, all 45 coefficient identities, and an explicit determinant law."""
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, combinations_with_replacement, product
from math import prod
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def local(name,path):
    spec = spec_from_file_location(name, ROOT/path)
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


double = local("double_pairs","double-paired-contingencies/check.py")
colors = local("three_colors","three-color-coupling/coupling.py")
Q,e = colors.Q,colors.e
P,PAIRS = double.PERMS,double.PAIRS
N = tuple(range(6))


def require(ok, msg):
    if not ok:
        raise ValueError(msg)


def cross(z,w):
    return z.real*w.imag-z.imag*w.real


def bary(p,vertices):
    if len(vertices)==1:
        require(all(x==vertices[0] for x in p),"not a point")
        return [[F(1)] for _ in p]
    if len(vertices)==2:
        x,y = vertices
        d = y-x
        require(d!=Q(0) and all(cross(d,z-x)==0 for z in p),"not on segment line")
        pars = [((z-x).conjugate()*d).real/colors.abs2(d) for z in p]
        weights = [[1-t,t] for t in pars]
    else:
        x,y,z = vertices
        den = cross(y-x,z-x)
        require(den!=0,"collinear triangle")
        weights = [[cross(y-t,z-t)/den, cross(z-t,x-t)/den, cross(x-t,y-t)/den] for t in p]
    require(all(min(row)>=0 and sum(row)==1 and sum(v*w for v,w in zip(vertices,row))==x
                for row,x in zip(weights,p)),"outside selected hull")
    return weights


def products(a,b):
    return [prod(a[i]+b[j] for i,j in enumerate(perm)) for perm in P]


def positive_cone():
    cells = tuple(product(range(3),range(3,6)))
    primitive_rows = {double.primitive(row) for row in double.monomial_rows().values()}
    generators = []
    for g,h in combinations_with_replacement(cells,2):
        vector = [int(i in g and j in h)+int(i in h and j in g) for i,j in PAIRS]
        generators.append(double.primitive(vector))
    require(len(set(generators))==45 and set(generators)<=primitive_rows,
            "45 bilinear forms must be positive multiples of Cauchy–Binet rows")


def main():
    positive_cone()
    a = (Q(-2,-2),Q(-1,3),Q(3,3))
    b = (Q(-1,-3),Q(-3,-1),Q(-3,-2))
    p = products(a,b)
    require(p==[Q(-14,22),Q(46,14),Q(-16,20),Q(-34,34),Q(-6,10),Q(56,12)],"products")
    vertices = p[3:]
    coords = bary(p,vertices)
    ref = [[95,80,18],[18,5,170],[80,108,5]]
    require(coords[:3]==[[F(x,193) for x in row] for row in ref], "inner barycentric rows")
    require(all(sum(col[:3])==sum(col[3:]) for col in zip(*coords)),"balanced affine coordinates")
    x,y,z = vertices
    den = cross(y-x,z-x)
    require(den==1544 and min(w*den for row in coords[:3] for w in row)==40,"strict area slack")
    require(max(colors.abs2(t+u) for t,u in product(a,b))<100 and
            max(colors.abs2(t-u) for t,u in combinations(p,2))<100**2,"distance bounds")
    eps = 6*12**2*F(1,16384)
    require(eps<=F(1,16) and 400*eps+4*eps*eps<40,"uniform perturbation proof")
    # Independent direct determinant on the repository's arbitrary dense complex unitary.
    u = colors.fixture()
    colors.check_unitary(u)
    require(all(z!=Q(0) for row in u for z in row),"dense test matrix")
    def table(c):
        return tuple(tuple(sum(c[i]==j for i in (2*r,2*r+1)) for j in range(3))
                     for r in range(3))
    ix = {table(tuple(val for r in range(3) for val in (P[i][r],P[j][r]))):k
          for k,(i,j) in enumerate(PAIRS)}
    q = [F(0)]*21
    for c,w in zip(colors.COLORINGS,colors.canon(u)):
        q[ix[table(c)]] += w
    weights = {}
    for k,l in combinations_with_replacement(range(3),2):
        value = sum(w*(coords[i][k]*coords[j][l]+coords[j][k]*coords[i][l])/2
                    for w,(i,j) in zip(q,PAIRS))
        weights[k,l] = value*(1 if k==l else 2)
    require(min(weights.values())>=0 and sum(weights.values())==1,"six actual vertex weights")
    db = e.diagonal([b[i//2] for i in N])
    m = e.multiply(e.multiply(u,db),e.adjoint(u))
    for i in N:
        m[i][i] += a[i//2]
    target = e.determinant(m)
    require(target.imag!=0 and target==sum(w*vertices[k]*vertices[l] for (k,l),w in weights.items()),
            "original determinant recovered")
    require(bary([Q(2)]*6,[Q(2)])==[[1]]*6,"singleton")
    require(bary([Q(0),Q(1),Q(2)],[Q(0),Q(2)])[1]==[F(1,2),F(1,2)],"segment")
    print("All 45 positive forms and the three barycentric dimensions pass.")
    print("The complex triangle has area determinant 1544 and inner margin 40; radius 1/16384.")
    print("The six nonnegative product weights recover the dense complex unitary determinant.")


if __name__ == "__main__":
    main()
