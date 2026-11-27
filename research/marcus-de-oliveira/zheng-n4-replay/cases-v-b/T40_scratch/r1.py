import sys, time, itertools
sys.path.insert(0,'C:/Users/nrvig/Projects/math-commons/research/marcus-de-oliveira/zheng-n4-replay/cases-v-b')
import replay_T38_common as c
from flint import fmpz_mpoly_ctx, Ordering, fmpz
from fractions import Fraction
names=['q','p','j','c','e','d','b','a','z']
ctx=fmpz_mpoly_ctx.get(names, Ordering.lex)
q,p,j,cc,e,d,b,a,z=ctx.gens()
src={5:q,1:p,12:j,9:cc,2:e,35:d,33:b,29:a}
u=[ctx.constant(0)]*14
for f,v in src.items():
    for k in range(14):
        x=c.F_COORD[f][k]
        assert x.denominator==1
        if x!=0: u[k]=u[k]+int(x)*v
u[13]=u[13]+z
print([str(x) for x in u])
C=c.cubic_C(u)
print(C.to_dict() and len(C.terms()))
# C0, Q
C0=C.subs({'z':0}) if False else None
