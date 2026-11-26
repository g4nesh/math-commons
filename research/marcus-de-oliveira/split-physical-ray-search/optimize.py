"""Bounded approximate search on actual unitary matrices, not relaxed tables.

Requires NumPy and SciPy for discovery; values are not positivity proofs.
"""
from pathlib import Path
from itertools import permutations
import argparse
import json
import random
import time
import numpy as np
import scipy
from scipy.optimize import minimize

HERE=Path(__file__).resolve().parent
ROW=(0,0,1,1,2,3)
COL=(0,0,1,1,2,2)
PERMS=np.array(list(permutations(range(6))))
SIGNS=np.array([(-1)**sum(p[i]>p[j] for i in range(6) for j in range(i+1,6)) for p in PERMS])
TABLES=sorted({tuple(sum(ROW[i]==r and COL[p[i]]==s for i in range(6)) for r in range(4) for s in range(3)) for p in PERMS})
LOOKUP={t:i for i,t in enumerate(TABLES)}
INDICES=np.array([LOOKUP[tuple(sum(ROW[i]==r and COL[p[i]]==s for i in range(6)) for r in range(4) for s in range(3))] for p in PERMS])
OFF=[(i,j) for i in range(6) for j in range(i+1,6)]


def hermitian(theta):
    h=np.diag(np.asarray(theta[:6],dtype=complex))
    for k,(i,j) in enumerate(OFF):
        h[i,j]=theta[6+2*k]+1j*theta[7+2*k]
        h[j,i]=h[i,j].conjugate()
    return h


def fixture(vectors):
    u=np.eye(6,dtype=complex)
    for record in vectors:
        v=np.array([complex(a,b) for a,b in record])
        u=(np.eye(6)-2*np.outer(v,v.conj())/(v.conj()@v))@u
    return u


def unitary(theta,base):
    inv=np.linalg.inv(np.eye(6)+1j*hermitian(theta))
    return (2*inv-np.eye(6))@base,inv


def objective(theta,base,ray):
    u,inv=unitary(theta,base)
    entries=u[np.arange(6)[:,None],PERMS.T]
    multipliers=SIGNS*np.asarray(ray)[INDICES]/sum(ray)
    products=np.prod(entries,axis=0)
    numerator=np.dot(multipliers,products)
    det=np.linalg.det(u)
    deriv=np.zeros((6,6),dtype=complex)
    for i in range(6):
        reduced=np.prod(np.delete(entries,i,axis=0),axis=0)
        np.add.at(deriv[i],PERMS[:,i],multipliers*reduced)
    deriv=(deriv-numerator*np.linalg.inv(u).T)/det
    # dg=Tr(deriv^T dU), dU=-2i inv*dH*inv*base.
    m=-2j*inv@base@deriv.T@inv
    gradient=[m[i,i].real for i in range(6)]
    for i,j in OFF:
        gradient.extend(((m[j,i]+m[i,j]).real,(1j*(m[j,i]-m[i,j])).real))
    return (numerator/det).real,np.array(gradient)


def gradient_check(base,ray):
    theta=np.random.default_rng(177).normal(0,.1,36)
    value,analytic=objective(theta,base,ray)
    step=1e-6
    numerical=[]
    for j in range(36):
        d=np.zeros(36);d[j]=step
        numerical.append((objective(theta+d,base,ray)[0]-objective(theta-d,base,ray)[0])/(2*step))
    error=float(np.max(np.abs(analytic-numerical)))
    if not error<1e-8:
        raise ValueError('analytic gradient disagreement: '+str(error))
    return error


def run(output,rounds=4,maxiter=200,seconds=480):
    rays=json.loads((HERE/'rays.json').read_text())
    rng=random.Random(202611260634)
    vectors=[[[[rng.randrange(-3,4),rng.randrange(-3,4)] for _ in range(6)] for _ in range(4)] for _ in range(rounds)]
    bases=[fixture(x) for x in vectors]
    gradient_error=gradient_check(bases[0],rays['records'][0]['ray'])
    result={'schema':'mdo.split-physical-ray-search.approximate.v1','seed':202611260634,
            'budget':{'rounds':rounds,'classes':len(rays['records']),'maxiter':maxiter,'wall_seconds':seconds},
            'versions':{'numpy':np.__version__,'scipy':scipy.__version__},
            'gradient_max_abs_error':gradient_error,'start_vectors':vectors,'runs':[],
            'negative_threshold':-1e-8,'interpretation':'finite approximate search, not a universal sign certificate'}
    start=time.monotonic()
    for repeat,base in enumerate(bases):
        for record in rays['records']:
            if time.monotonic()-start>=seconds:
                result['stop_reason']='wall budget'
                output.write_text(json.dumps(result,indent=2)+'\n')
                return result
            ray=record['ray']
            initial,_=objective(np.zeros(36),base,ray)
            answer=minimize(objective,np.zeros(36),args=(base,ray),jac=True,method='L-BFGS-B',
                            options={'maxiter':maxiter,'gtol':1e-10,'ftol':1e-14,'maxls':40})
            value,gradient=objective(answer.x,base,ray)
            u,_=unitary(answer.x,base)
            row={'orbit':record['id'],'start':repeat,'initial_value':float(initial),'value':float(value),
                 'theta':answer.x.tolist(),'iterations':int(answer.nit),'evaluations':int(answer.nfev),
                 'success':bool(answer.success),'status':int(answer.status),'message':str(answer.message),
                 'gradient_inf':float(np.max(np.abs(gradient))),
                 'unitarity_max_error':float(np.max(np.abs(u.conj().T@u-np.eye(6)))),
                 'elapsed_seconds':time.monotonic()-start}
            result['runs'].append(row)
            result['stop_reason']='negative candidate' if value<-1e-8 else 'running'
            output.write_text(json.dumps(result,indent=2)+'\n')
            print(repeat,record['id'],format(value,'.9g'),answer.nit,answer.success,flush=True)
            if value<-1e-8:
                return result
    result['stop_reason']='all scheduled runs completed'
    output.write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=HERE/'search.json')
    parser.add_argument('--rounds',type=int,default=4)
    parser.add_argument('--maxiter',type=int,default=200)
    parser.add_argument('--seconds',type=float,default=480)
    args=parser.parse_args()
    if not 1<=args.rounds<=4 or not 1<=args.maxiter<=200 or not 0<args.seconds<=480:
        parser.error('bounded budget exceeded')
    run(args.output,args.rounds,args.maxiter,args.seconds)
