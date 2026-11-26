"""Exact complex matrix checks in Q(i,zeta_3), without floating point."""
from fractions import Fraction as F
from functools import cache
from itertools import combinations
from model import SITES,ROWS,coefficients,require

ZERO=(F(0),)*4
ONE=(F(1),F(0),F(0),F(0))
I=(F(0),F(0),F(1),F(0))
ZETA=(F(0),F(1),F(0),F(0))
# sqrt(3) = -i*(1+2*zeta)
ROOT3=(F(0),F(0),F(-1),F(-2))

def add(z,w):return tuple(a+b for a,b in zip(z,w))
def scale(z,x):return tuple(a*x for a in z)
def emul(z,w):
    a,b=z;c,d=w
    return a*c-b*d,a*d+b*c-b*d

def mul(z,w):
    ac=emul(z[:2],w[:2]);bd=emul(z[2:],w[2:])
    ad=emul(z[:2],w[2:]);bc=emul(z[2:],w[:2])
    return ac[0]-bd[0],ac[1]-bd[1],ad[0]+bc[0],ad[1]+bc[1]

def conj(z):
    a,b,c,d=z
    return a-b,-b,-c+d,d

def total(values):
    answer=ZERO
    for v in values:answer=add(answer,v)
    return answer

def gaussian(a,b=0):return F(a),F(0),F(b),F(0)

def determinant3(rows):
    a,b,c=rows
    return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])

def replay_fixture(ab):
    blocks=[];axes=[]
    for a,b in ab:
        require(add(mul(a,conj(a)),mul(b,conj(b)))==ONE,'unit quaternion row')
        blocks.append(((a,b),(scale(conj(b),-1),conj(a))))
        cross=mul(conj(a),b)
        require(cross[1]==cross[3]==0,'Gaussian axis coefficients')
        difference=add(mul(a,conj(a)),scale(mul(b,conj(b)),-1))
        axes.append((2*cross[0],-2*cross[2],difference[0]))
    roots=(ONE,ZETA,mul(ZETA,ZETA))
    M=[[mul(blocks[i//2][i%2][j%2],roots[(i//2)*(j//2)%3]) for j in SITES] for i in SITES]
    require(all(total(mul(M[i][k],conj(M[j][k])) for k in SITES)==scale(ONE,3 if i==j else 0)
                for i in SITES for j in SITES),'exact unitarity with omitted sqrt(3)')
    @cache
    def det(rows,cols):
        if not rows:return ONE
        return total(scale(mul(M[rows[0]][c],det(rows[1:],cols[:j]+cols[j+1:])),(-1)**j)
                     for j,c in enumerate(cols))
    dots=[sum(a*b for a,b in zip(axes[i],axes[j])) for i,j in ((0,1),(0,2),(1,2))]
    tau=determinant3(axes)
    coordinates=[ONE]+[scale(ONE,x) for x in dots]+[scale(ROOT3,-tau)]
    for rows,h in ROWS:
        actual=ZERO
        for cols in combinations(SITES,len(rows)):
            if tuple(sum(c//2==j for c in cols) for j in range(3))==h:
                value=det(rows,cols)
                actual=add(actual,scale(mul(value,conj(value)),F(1,3**len(rows))))
        predicted=total(scale(x,q) for x,q in zip(coordinates,coefficients(rows,h)))
        require(actual==predicted,'raw squared-minor histogram does not match analytic formula')
    return {'pair_dots':[str(x) for x in dots],'triple_product':str(tau)}


def replay():
    fixtures=[[(gaussian(1),ZERO)]*3,
              [(gaussian(1),ZERO),(gaussian(F(3,5)),gaussian(F(4,5))),
               (gaussian(F(5,13)),gaussian(0,F(12,13)))],
              [(gaussian(F(1,3),F(2,3)),gaussian(F(2,3))),
               (gaussian(F(3,5)),gaussian(0,F(4,5))),
               (gaussian(F(5,13)),gaussian(F(12,13)))]]
    return [replay_fixture(f) for f in fixtures]

if __name__=='__main__':
    result=replay()
    print('PASS: 3 exact unitary fixtures, 924 minors and 358 grouped coefficients each.')
    print(result)
