"""Historical exact slice grid that missed the analytic negative-g0 family.
Finite misses provide no evidence of a universal inequality.
"""
import importlib.util,itertools,json
from pathlib import Path
from fractions import Fraction as F
p=Path(__file__).with_name('verify.py')
s=importlib.util.spec_from_file_location('b',p);b=importlib.util.module_from_spec(s);s.loader.exec_module(b)
pairs=[(F(3,5),F(4,5)),(F(5,13),F(12,13))]
mixes=[(F(0),F(1)),(F(3,5),F(4,5)),(F(1),F(0))]
phases=[b.z(1),b.z(0,1),b.z(F(3,5),F(4,5))]
for kind in ['x1zero','x2zero','x4zero']:
 minimum=None;minA=None;viol=0;worst=None;count=0
 for ip,ia,ib,ir,phase in itertools.product(range(2),range(3),range(3),range(3),range(3)):
  p,q=pairs[ip];C,S=(F(0),F(1)) if kind=='x1zero' else ((F(1),F(0)) if kind=='x2zero' else pairs[0]);E,T=(F(1),F(0)) if kind=='x4zero' else pairs[0]
  ra,sa=mixes[ia];rb,sb=mixes[ib];ar,rr=mixes[ir]
  U=b.boundary(p,q,C,S,E,T,b.scale(ra,phases[phase]),b.scale(sa,phases[(phase+1)%3]),b.scale(rb,phases[(phase+2)%3]),b.z(sb),b.z(ar),rr)
  assert b.unitary(U);A,B,g=b.quantities(U);m=g+A
  violated=m<0 and m*m>4*A*B
  viol+=violated;count+=1
  if minimum is None or g<minimum:minimum=g;worst=b.fixture(U,kind)
  if minA is None or m<minA:minA=m
 print(json.dumps({'slice':kind,'count':count,'minimum_g0':str(minimum),'minimum_g0_plus_A':str(minA),'star_violations':viol}))
 if minimum<0:print(json.dumps(worst))
