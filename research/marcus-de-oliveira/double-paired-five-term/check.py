"""Exact five-term reductions and obstruction to the Hodge-only relaxation."""
from fractions import Fraction as F
from itertools import combinations,permutations,product
from functools import cache
import json

P3=((0,1,2),(1,2,0),(2,0,1),(0,2,1),(1,0,2),(2,1,0))
PAIRS=tuple(combinations(range(6),2))
QPAIRS=tuple((i,j) for i in range(6) for j in range(i,6))
S3=tuple(combinations(range(6),3))
TARGET={(1,2),(1,5),(2,4),(4,4),(5,5)}

def require(ok,message):
    if not ok:raise ValueError(message)

def sign(p):return (-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))
def hist(s):return tuple(sum(i//2==j for i in s) for j in range(3))
def complement(s):return tuple(i for i in range(6) if i not in s)
def hodge_sign(s):return sign(s+complement(s))

def profile_identity():
    rows={}
    for k,(i,j) in enumerate(QPAIRS):
        assigned=P3[i]+P3[j]
        for take in product((False,True),repeat=6):
            a=tuple(sum(not take[v] and v%3==r for v in range(6)) for r in range(3))
            b=tuple(sum(take[v] and assigned[v]==r for v in range(6)) for r in range(3))
            rows.setdefault((a,b),[0]*21)[k]+=1
    terms=[((1,2,0),(1,1,1),-1),((1,1,1),(1,2,0),1),
           ((2,1,0),(2,1,0),2),((2,0,1),(2,0,1),2)]
    actual=[sum(c*rows[(a,b)][j] for a,b,c in terms) for j in range(21)]
    require(actual==[4*int(p in TARGET) for p in QPAIRS],'21-coordinate third-compound identity')

# Sparse commutative polynomials: monomial is a sorted tuple of variable ids.
def plus(a,b):
    c=dict(a)
    for m,v in b.items():c[m]=c.get(m,0)+v
    return {m:v for m,v in c.items() if v}
def times(a,b):
    c={}
    for m,x in a.items():
        for n,y in b.items():
            k=tuple(sorted(m+n));c[k]=c.get(k,0)+x*y
    return {m:v for m,v in c.items() if v}
def scaled(a,c):return {m:c*x for m,x in a.items() if c*x}
def detp(M):
    answer={}
    for p in permutations(range(len(M))):
        term={():sign(p)}
        for i,j in enumerate(p):term=times(term,M[i][j])
        answer=plus(answer,term)
    return answer
def mp(A,B):return [[sum_poly(times(x,y) for x,y in zip(row,col)) for col in zip(*B)] for row in A]
def sum_poly(polys):
    out={}
    for p in polys:out=plus(out,p)
    return out
def adj(A):return [[A[1][1],scaled(A[0][1],-1)],[scaled(A[1][0],-1),A[0][0]]]

def block_identity():
    blocks=[[[{(4*k+2*i+j,):1} for j in range(2)] for i in range(2)] for k in range(6)]
    A,B,C,D,E,G=blocks
    zero=[[{},{}],[{},{}]]
    grid=[[zero,A,B],[C,zero,D],[E,G,zero]]
    M=[[grid[i//2][j//2][i%2][j%2] for j in range(6)] for i in range(6)]
    cycle=mp(mp(mp(mp(mp(adj(B),A),adj(G)),E),adj(C)),D)
    trace=plus(cycle[0][0],cycle[1][1])
    diagonal=plus(times(times(detp(A),detp(D)),detp(E)),times(times(detp(B),detp(C)),detp(G)))
    require(detp(M)==plus(diagonal,trace),'universal six-block trace identity')
    require(detp(M)!=plus(diagonal,scaled(trace,-1)),'changed cycle sign must fail')
    return len(trace),len(detp(M))

# Literal relaxed map on the 20 triples; each complement pair is mapped together.
IMAGE=((0,3,5),(0,3,4),(0,1,4),(0,1,5),(0,1,2),(0,2,3),(0,2,4),(0,2,5),(0,4,5),(0,1,3),
       (2,4,5),(1,2,3),(1,3,4),(1,3,5),(1,4,5),(3,4,5),(2,3,4),(2,3,5),(1,2,5),(1,2,4))

def relaxed_coefficient(r,c):
    a,b=hist(r),hist(c)
    return -int(a==(1,0,2) and b==(1,1,1))+int(a==(1,1,1) and b==(1,2,0))+2*int(a==(0,1,2) and b==(2,1,0))+2*int(a==(0,2,1) and b==(2,0,1))

def hodge_obstruction():
    mapping=dict(zip(S3,IMAGE))
    require(set(mapping.values())==set(S3),'basis permutation')
    require(all(mapping[complement(s)]==complement(mapping[s]) for s in S3),'complement pairing')
    signs={}
    for s in S3:
        if 0 in s:
            signs[s]=1
            signs[complement(s)]=hodge_sign(mapping[s])*hodge_sign(s)
    C=[[int(r==mapping[c])*signs[c] for c in S3] for r in S3]
    J=[[int(r==complement(c))*hodge_sign(c) for c in S3] for r in S3]
    mul=lambda A,B:[[sum(x*y for x,y in zip(row,col)) for col in zip(*B)] for row in A]
    require(mul(C,list(zip(*C)))==[[int(i==j) for j in range(20)] for i in range(20)],'relaxed unitarity')
    require(mul(C,J)==mul(J,C),'Hodge commutation')
    require(mul(mul(C,J),list(zip(*C)))==J,'symplectic form')
    value=sum(relaxed_coefficient(r,c)*C[i][j]**2 for i,r in enumerate(S3) for j,c in enumerate(S3))
    require(value==-2,'relaxed functional')
    image={mapping[s]:signs[s] for s in [(0,1,2),(0,1,4)]}
    p=lambda s:image.get(s,0)
    plucker=p((0,1,3))*p((0,4,5))-p((0,1,4))*p((0,3,5))+p((0,1,5))*p((0,3,4))
    require(plucker==-1,'nondecomposable image')
    return {'four_phi':value,'phi':'-1/2','Plucker_residual':plucker,'is_actual_U6_counterexample':False}

def replay():
    profile_identity()
    counts=block_identity()
    result=hodge_obstruction()
    return {'status':'PASS','profile_coordinates':21,'block_trace_monomials':counts[0],
            'off_diagonal_determinant_monomials':counts[1],'relaxation':result,
            'universal_five_term_sign':'OPEN'}


# Exact physical fixtures check index conventions; they are not sign proofs.
class Gaussian:
    def __init__(self,a=0,b=0):self.a,self.b=F(a),F(b)
    def __add__(self,z):
        if not isinstance(z,Gaussian):z=Gaussian(z)
        return Gaussian(self.a+z.a,self.b+z.b)
    __radd__=__add__
    def __neg__(self):return Gaussian(-self.a,-self.b)
    def __sub__(self,z):return self+-z
    def __mul__(self,z):
        if not isinstance(z,Gaussian):z=Gaussian(z)
        return Gaussian(self.a*z.a-self.b*z.b,self.a*z.b+self.b*z.a)
    __rmul__=__mul__
    def __eq__(self,z):
        if not isinstance(z,Gaussian):z=Gaussian(z)
        return (self.a,self.b)==(z.a,z.b)
    def conjugate(self):return Gaussian(self.a,-self.b)
    def norm(self):return self.a*self.a+self.b*self.b


def physical_fixtures():
    tables={tuple(sum(P3[i][r]==c for i in pair) for r in range(3) for c in range(3)):pair for pair in QPAIRS}
    require(len(tables)==21,'contingency decomposition')
    results=[]
    for seed in range(2):
        U=[[Gaussian(int(i==j)) for j in range(6)] for i in range(6)]
        for k,(i,j) in enumerate([(0,2),(1,5),(2,4),(0,3),(1,4),(3,5),(0,1),(2,5)]):
            t=F(k+seed+1,k+seed+4);c=(1-t*t)/(1+t*t);s=2*t/(1+t*t)
            phase=Gaussian(F(3,5),F(4,5)) if k%2 else Gaussian(F(5,13),F(12,13))
            a,b=U[i][:],U[j][:]
            U[i]=[c*x+s*phase*y for x,y in zip(a,b)]
            U[j]=[-s*phase.conjugate()*x+c*y for x,y in zip(a,b)]
        require(all(sum(U[i][k]*U[j][k].conjugate() for k in range(6))==int(i==j)
                    for i in range(6) for j in range(6)),'physical unitary fixture')
        @cache
        def minor(rows,cols):
            if not rows:return Gaussian(1)
            return sum((-1)**j*U[rows[0]][c]*minor(rows[1:],cols[:j]+cols[j+1:]) for j,c in enumerate(cols))
        z=minor(tuple(range(6)),tuple(range(6))).conjugate()
        q={pair:F(0) for pair in QPAIRS}
        for p in permutations(range(6)):
            term=Gaussian(sign(p))
            for i,j in enumerate(p):term=term*U[i][j]
            table=tuple(sum(i//2==r and p[i]//2==c for i in range(6)) for r in range(3) for c in range(3))
            q[tables[table]]+=(z*term).a
        require(sum(q.values())==1,'canonical law normalization')
        y=sum(q[p] for p in TARGET)
        four_phi=sum(relaxed_coefficient(r,c)*minor(r,c).norm() for r in S3 for c in S3)
        require(4*y==four_phi,'physical third-compound bridge')
        blocks={(r,c):minor((2*r,2*r+1),(2*c,2*c+1)) for r in range(3) for c in range(3)}
        require(y==blocks[2,0].norm()+blocks[1,0].norm()+q[1,2]-q[1,1]-q[2,2],
                'physical two-block-norm bridge')
        for i in [1,2]:
            value=Gaussian(1)
            for r in range(3):value=value*blocks[r,P3[i][r]]
            require(q[i,i]==(z*value).a,'diagonal block coefficient')
        results.append({'seed':seed,'y':str(y),'negative_starting_coefficients':sum(v<0 for v in q.values())})
    return results

if __name__=='__main__':
    report=replay()
    report['physical_fixtures']=physical_fixtures()
    print(json.dumps(report,indent=2))
