import pickle, time
from flint import fmpz_mpoly_ctx, Ordering
names=['q','p','j','c','e','d','b','a','z']
ctx=fmpz_mpoly_ctx.get(names, Ordering.lex)
Rs={k:ctx.from_dict(v) for k,v in pickle.load(open('Rs.pkl','rb')).items()}
# wall j=1: p=P/(1+P), multiply by (1+P)^2
ctx3=fmpz_mpoly_ctx.get(['q','P','e','b','a'],Ordering.lex)
q3,P3,e3,b3,a3=ctx3.gens()
for nm,R in Rs.items():
    degs=R.degrees(); dp=degs[1]
    out=ctx3.constant(0)
    for ex,co in R.to_dict().items():
        kq,kp,kj,_,ke,_,kb,ka,_=ex
        out=out+int(co)*q3**kq*e3**ke*b3**kb*a3**ka*P3**kp*(1+P3)**(dp-kp)
    # j=1: substitute j=1 (j^kj -> 1) done implicitly (kj ignored)
    dct=out.to_dict(); cs=list(dct.values())
    print('wall',nm,'degp',dp,len(cs),all(x<0 for x in cs),dct.get((0,)*5))
# high-p: p=1+P, j=J/((1+P)(1+J)); multiply by [(1+P)(1+J)]^degj
ctx4=fmpz_mpoly_ctx.get(['q','P','J','e','b','a'],Ordering.lex)
q4,P4,J4,e4,b4,a4=ctx4.gens()
def hp(R):
    dj=R.degrees()[2]
    out=ctx4.constant(0)
    for ex,co in R.to_dict().items():
        kq,kp,kj,_,ke,_,kb,ka,_=ex
        out=out+int(co)*q4**kq*e4**ke*b4**kb*a4**ka*(1+P4)**kp*J4**kj*((1+P4)*(1+J4))**(dj-kj)
    return out
for nm in ['gamma','ABb-Nja']:
    S=hp(Rs[nm]); dct=S.to_dict(); cs=list(dct.values())
    print('highp',nm,len(cs),all(x<0 for x in cs),dct.get((0,)*6))
# sum certificate AB(F+ + F-) + 2 N alpha = R3 + R4 - ... R3+R4 = AB(F++F-)+2Nal
S=hp(Rs['ABF++Na']+Rs['ABF-+Na']); dct=S.to_dict(); cs=list(dct.values())
print('highp sum',len(cs),all(x<0 for x in cs),dct.get((0,)*6))
# the individual ones in high-p: need not be negative
for nm in ['ABF++Na','ABF-+Na']:
    S=hp(Rs[nm]); cs=list(S.to_dict().values()); print('highp indiv',nm,len(cs),sum(1 for x in cs if x>0),'positive coeffs')
