from fractions import Fraction as F
from itertools import product
import argparse, json, platform, time
from pathlib import Path

def z(x=0,y=0): return F(x),F(y)
def add(a,b): return a[0]+b[0],a[1]+b[1]
def neg(a): return -a[0],-a[1]
def sub(a,b): return add(a,neg(b))
def conj(a): return a[0],-a[1]
def mul(a,b): return a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]
def scale(t,a): return t*a[0],t*a[1]
def norm(a): return a[0]**2+a[1]**2
def vs(a,b): return [add(x,y) for x,y in zip(a,b)]
def vm(t,a): return [mul(t,x) for x in a]
def dot(a,b): return sum((mul(x,conj(y))[0] for x,y in zip(a,b)),F(0)),sum((mul(x,conj(y))[1] for x,y in zip(a,b)),F(0))
def minor(U,rs,cs):
 i,j=rs;k,l=cs
 return sub(mul(U[i][k],U[j][l]),mul(U[i][l],U[j][k]))
def D(U,rs,cs): return norm(minor(U,rs,cs))
def B(U,i,j): return norm(U[i][j])
def quantities(U):
 a=D(U,(0,1),(0,1)); b=D(U,(0,2),(0,3))
 g=B(U,0,0)+B(U,2,2)+D(U,(0,3),(2,3))-D(U,(0,1),(0,3))
 assert g==D(U,(0,2),(0,2))+D(U,(0,3),(0,3))+D(U,(0,2),(1,2))+b-D(U,(0,1),(2,3))
 return a,b,g

def boundary(p,q,C,S,E,T,beta,gamma,delta,eta,alpha=None,r=None):
 x=[z(p*C),z(p*S),z(q*E),z(q*T)]
 f=[z(S),z(-C),z(),z()]
 g=[z(),z(),z(T),z(-E)]
 h=[z(q*C),z(q*S),z(-p*E),z(-p*T)]
 a=vs(vm(beta,g),vm(gamma,h))
 v=vs(vm(neg(conj(gamma)),g),vm(conj(beta),h))
 if alpha is not None:
  f,a=vs(vm(z(r),f),vm(neg(conj(alpha)),a)),vs(vm(alpha,f),vm(z(r),a))
 b=vs(vm(delta,f),vm(eta,v))
 c=vs(vm(neg(conj(eta)),f),vm(conj(delta),v))
 return [x,a,b,c]

def unitary(U):
 return all(dot(U[i],U[j])==z(i==j) for i,j in product(range(4),repeat=2))

def fixture(U,kind):
 a,b,g=quantities(U);P=B(U,0,0)+B(U,0,1);K=B(U,0,0)*B(U,0,1)*B(U,0,3)
 return {'kind':kind,'schema':'exact-gaussian-rational-matrix-v1',
         'U':[[[str(re),str(im)] for re,im in row] for row in U],
         'A':str(a),'C':str(b),'g0':str(g),'P':str(P),'K':str(K),
         'condition_slack_K2_minus_49PA':str(K*K-49*P*a),
         't_star_squared':str(a/b) if b else None}

def run():
 start=time.process_time()
 pairs=[(F(3,5),F(4,5)),(F(5,13),F(12,13)),(F(8,17),F(15,17))]
 angles=[(F(0),F(1)),(F(3,5),F(4,5)),(F(1),F(0))]
 mixes=[(F(0),F(1)),(F(3,5),F(4,5)),(F(5,13),F(12,13)),(F(1),F(0))]
 phases=[z(1),z(0,1),z(F(3,5),F(4,5)),z(F(-5,13),F(12,13))]
 count=complex_count=zero_gap=0; min_gap=None; fixtures=[]
 for ip,it,jt,ia,ib in product(range(3),range(3),range(3),range(4),range(4)):
  p,q=pairs[ip];C,S=angles[it];E,T=angles[jt];ra,sa=mixes[ia];rb,sb=mixes[ib]
  beta=scale(ra,phases[(ia+ip)%4]);gamma=scale(sa,phases[(ia+jt+1)%4])
  delta=scale(rb,phases[(ib+it+1)%4]);eta=scale(sb,phases[(ib+2)%4])
  U=boundary(p,q,C,S,E,T,beta,gamma,delta,eta)
  assert unitary(U)
  a,b,g0=quantities(U);assert a==0
  u,v,w,r=p*C,p*S,q*E,q*T
  A1=(u*u/(q*q))*norm(sub(scale(r,beta),scale(w/p,gamma)))+(q*q*v*v/(p*p))*norm(gamma)
  A2=u*u+(v*v/(q*q))*norm(add(scale(w,beta),scale(r/p,gamma)))
  assert g0==norm(delta)*A1+norm(eta)*A2
  gap=u*u*v*v*r*r/(p*p)
  assert g0>=gap>=0
  assert A1>=gap and A2>=gap
  count+=1;complex_count+=any(im for row in U for re,im in row);zero_gap+=g0==0
  if min_gap is None or g0-gap<min_gap:min_gap=g0-gap
  if not fixtures and ip==0 and it==1 and jt==1 and ia==1 and ib==1:
   fixtures.append(fixture(U,'nondegenerate exact boundary'))

 # Perturbations use exact unitary rotations away from A=0. No square roots or
 # rounded input are used in checking the polynomial sufficient condition.
 perturbations=positive_A=sufficient=uncovered_ratio=0
 for ip,ia,ib,k in product(range(3),range(4),range(4),[F(1,1000),F(1,100000),F(1,5)]):
  p,q=pairs[ip];C,S=angles[1];E,T=angles[1];ra,sa=mixes[ia];rb,sb=mixes[ib]
  beta=scale(ra,phases[(ia+ip)%4]);gamma=scale(sa,phases[(ia+2)%4])
  delta=scale(rb,phases[(ib+2)%4]);eta=scale(sb,phases[(ib+1)%4])
  alpha=scale(2*k/(1+k*k),phases[(ia+ib)%4]);r=(1-k*k)/(1+k*k)
  U=boundary(p,q,C,S,E,T,beta,gamma,delta,eta,alpha,r)
  assert unitary(U)
  A,Cminor,g0=quantities(U);P=p*p;K=B(U,0,0)*B(U,0,1)*B(U,0,3);gap=K/P
  assert A==P*norm(alpha) and 0<A<P
  if g0<gap:assert (gap-g0)**2<=49*A/P
  positive_A+=A>0;perturbations+=1
  if 49*P*A<=K*K:
   assert g0>=0;sufficient+=1
   if Cminor and 25*A<Cminor:
    uncovered_ratio+=1
    if len(fixtures)==1:fixtures.append(fixture(U,'nonzero small minor below t=1/5 with sufficient condition'))

 # p=0: use a general complex three-row frame, not the nondegenerate formula.
 pzero=qzero=0
 for i,j in product(range(4),repeat=2):
  w,zz=F(3,5),F(4,5)
  H=[z(),z(),z(zz),z(-w)]
  rows=[[z(1),z(),z(),z()],[z(),z(1),z(),z()],H]
  for turn,(left,right) in enumerate([(0,1),(1,2),(0,2)]):
   cr,sr=mixes[(i+turn)%4];ss=scale(sr,phases[(j+turn)%4])
   L,R=rows[left],rows[right]
   rows[left]=vs(vm(z(cr),L),vm(ss,R))
   rows[right]=vs(vm(neg(conj(ss)),L),vm(z(cr),R))
  U=[[z(),z(),z(w),z(zz)]]+rows
  assert unitary(U)
  A,_,g0=quantities(U)
  assert A==0 and g0==zz*zz*norm(U[1][1])+w*w*norm(dot(U[3],H)) and g0>=0
  pzero+=1
  beta=scale(mixes[i][0],phases[i]);gamma=scale(mixes[i][1],phases[(i+1)%4])
  delta=scale(mixes[j][0],phases[j]);eta=scale(mixes[j][1],phases[(j+2)%4])
  U=boundary(F(1),F(0),F(3,5),F(4,5),F(3,5),F(4,5),beta,gamma,delta,eta)
  assert unitary(U)
  A,_,g0=quantities(U)
  assert A==0 and g0==F(9,25)*norm(U[1][2])+norm(U[2][2]) and g0>=0
  qzero+=1

 report={'boundary_cases':count,'genuinely_complex_boundary_cases':complex_count,
  'boundary_zero_g0_cases':zero_gap,'min_boundary_g0_minus_gap':str(min_gap),
  'perturbation_cases':perturbations,'positive_A_perturbations':positive_A,
  'perturbations_satisfying_polynomial_condition':sufficient,
  'condition_cases_with_positive_C_and_t_star_below_1_over_5':uncovered_ratio,
  'p_zero_cases':pzero,'q_zero_cases':qzero,'total_exact_unitary_cases':count+perturbations+pzero+qzero,
  'failures':0,'python':platform.python_version(),'cpu_seconds':time.process_time()-start,
  'nature':'fixed exact Gaussian-rational regression; analytic proof is README sections2-4',
  'production_search':False}
 return report,fixtures

if __name__=='__main__':
 parser=argparse.ArgumentParser(description='Replay exact boundary and small-minor fixtures.')
 parser.add_argument('--write-reports',action='store_true',help='Rewrite report.json and fixtures.json beside this script.')
 args=parser.parse_args();report,fixtures=run()
 if args.write_reports:
  directory=Path(__file__).resolve().parent
  (directory/'report.json').write_text(json.dumps(report,indent=2)+'\n')
  (directory/'fixtures.json').write_text(json.dumps(fixtures,indent=2)+'\n')
 print(json.dumps(report,indent=2))
