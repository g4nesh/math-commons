import sys, time, itertools
sys.path.insert(0,'C:/Users/nrvig/Projects/math-commons/research/marcus-de-oliveira/zheng-n4-replay/cases-v-b')
import replay_T38_common as c
import sympy

def hessian_matrix(u, ctx):
    """14x14 Hessian of C at the point u (list of mpolys), entries linear."""
    ms, C, grad, hess = c.cubic_data()
    H=[[None]*14 for _ in range(14)]
    for i in range(14):
        for k in range(14):
            e=hess[i][k]
            tot=ctx.constant(0)
            if e!=0:
                for mono,coef in sympy.Poly(e,*ms).terms():
                    idx=[t for t,ee in enumerate(mono) if ee]
                    assert len(idx)==1 and mono[idx[0]]==1
                    tot=tot+int(coef)*u[idx[0]]
            H[i][k]=tot
    return H

def det_dp(M, rows=None, cols=None):
    """Exact determinant of a square matrix of fmpz_mpoly (Laplace DP over column subsets)."""
    n=len(M)
    ctxz=None
    states={0: None}
    first=True
    # choose to process rows in given order
    cur={0:None}
    for ri in range(n):
        nxt={}
        row=M[ri]
        nz=[k for k in range(n) if not row[k].is_zero()]
        for mask,val in cur.items():
            for k in nz:
                if mask>>k &1: continue
                inv=bin(mask>>(k+1)).count('1')
                term=row[k] if val is None else val*row[k]
                if inv&1: term=-term
                nm=mask|(1<<k)
                if nm in nxt: nxt[nm]=nxt[nm]+term
                else: nxt[nm]=term
        cur={m:v for m,v in nxt.items() if not v.is_zero()}
        if not cur: return None
    return cur[(1<<n)-1]
