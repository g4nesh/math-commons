import sys, pickle
sys.path.insert(0, r"C:\Users\nrvig\Projects\math-commons\research\marcus-de-oliveira\zheng-n4-replay\cases-v-b")
sys.dont_write_bytecode = True
import replay_T38_common as M
import flint
from flint import fmpz_mpoly_ctx, fmpz_mpoly
def ctxs(names):
    ctx = fmpz_mpoly_ctx.get(tuple(names), 'lex')
    return ctx, dict(zip(names, ctx.gens()))
def fromR(R, p, ctx):
    d = {}
    for k, c in p.items():
        d[tuple(R.exps(k))] = c
    return ctx.from_dict(d)
