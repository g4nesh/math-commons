"""Stdlib exact reconstruction of35 selected physical ray search endpoints.

These are finite signs, not a universal nonnegativity certificate.
"""
from fractions import Fraction as F
from pathlib import Path
import importlib.util
import json
import math
import sys

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('physical_ray_base',HERE.parent/'split-one-eigenvalue-pair'/'check.py')
base=importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)
Q,exact,require=base.Q,base.exact,base.require
OFF=[(i,j)for i in range(6)for j in range(i+1,6)]


def inverse(a):
    n=len(a)
    x=[row[:]+[Q(int(i==j))for j in range(n)]for i,row in enumerate(a)]
    for j in range(n):
        pivot=next((i for i in range(j,n)if x[i][j]!=Q(0)),None)
        require(pivot is not None,'invertible Cayley denominator')
        x[j],x[pivot]=x[pivot],x[j]
        v=x[j][j];x[j]=[a/v for a in x[j]]
        for i in range(n):
            if i!=j:
                v=x[i][j];x[i]=[a-v*b for a,b in zip(x[i],x[j])]
    return[row[n:]for row in x]


def rational_unitary(vectors,theta):
    require(len(theta)==36 and all(type(x) is int for x in theta),'36 integer chart numerators')
    h=[[Q(0)for _ in range(6)]for _ in range(6)]
    for i in range(6):h[i][i]=Q(F(theta[i],10**6))
    for k,(i,j)in enumerate(OFF):
        h[i][j]=Q(F(theta[6+2*k],10**6),F(theta[7+2*k],10**6))
        h[j][i]=h[i][j].conjugate()
    den=[[Q(int(i==j))+Q(0,1)*h[i][j]for j in range(6)]for i in range(6)]
    inv=inverse(den)
    cayley=[[2*inv[i][j]-Q(int(i==j))for j in range(6)]for i in range(6)]
    u=exact.multiply(cayley,base.fixture({'vectors':vectors}))
    base.check_unitary(u)
    return u


def prepare():
    search=json.loads((HERE/'search.json').read_text())
    rays=json.loads((HERE/'rays.json').read_text())
    result={'schema':'mdo.split-physical-ray-search.exact-endpoints.v1','chart_denominator':10**6,'records':[]}
    for ray in rays['records']:
        selected=min((r for r in search['runs']if r['orbit']==ray['id']),key=lambda r:r['value'])
        theta=[round(x*10**6)for x in selected['theta']]
        u=rational_unitary(search['start_vectors'][selected['start']],theta)
        q=base.canonical(u)
        value=sum(z*x for z,x in zip(ray['ray'],q))/sum(ray['ray'])
        result['records'].append({'orbit':ray['id'],'start':selected['start'],'theta_numerators':theta,'exact_value':str(value)})
        print('exact',ray['id'],'sign',(value>0)-(value<0),'approx',float(value),flush=True)
    (HERE/'exact-endpoints.json').write_text(json.dumps(result,indent=2)+'\n')


def verify():
    search=json.loads((HERE/'search.json').read_text())
    rays=json.loads((HERE/'rays.json').read_text())
    endpoints=json.loads((HERE/'exact-endpoints.json').read_text())
    require(search['schema']=='mdo.split-physical-ray-search.approximate.v1','search schema')
    require(endpoints['schema']=='mdo.split-physical-ray-search.exact-endpoints.v1'and endpoints['chart_denominator']==10**6,'exact endpoint schema')
    ids=rays['ids'];runs=search['runs']
    require(len(ids)==35 and len(set(ids))==35,'35 distinct orbit classes')
    require([(r['orbit'],r['start'])for r in runs]==[(i,s)for s in range(4)for i in ids],'all140 scheduled runs')
    require(all(r['iterations']<=200 and math.isfinite(r['value'])and all(math.isfinite(x)for x in r['theta'])for r in runs),'bounded finite records')
    require(all(r['value']>=search['negative_threshold']for r in runs),'no robust negative in archived search')
    require(search['gradient_max_abs_error']<1e-8,'archived gradient check')
    require([r['orbit']for r in endpoints['records']]==ids,'one exact endpoint per class')
    s=base.structure();values=[]
    for record,ray in zip(endpoints['records'],rays['records']):
        z=ray['ray']
        require(len(z)==33 and all(type(x)is int and x>=0 for x in z)and sum(z)>0,'positive integral ray')
        # Membership in the coefficient row span is independently checked by rank.
        require(len(base.independent_rows(s['matrix']+[z]))==21,'dual ray belongs to row space')
        selected=min((r for r in runs if r['orbit']==record['orbit']),key=lambda r:r['value'])
        require(record['start']==selected['start'] and record['theta_numerators']==[round(x*10**6)for x in selected['theta']],'declared selected quantization')
        u=rational_unitary(search['start_vectors'][record['start']],record['theta_numerators'])
        value=sum(x*y for x,y in zip(z,base.canonical(u)))/sum(z)
        require(value==F(record['exact_value']),'exact physical ray evaluation')
        require(value>=0,'negative exact endpoint found')
        values.append(value)
    return {'scope':'finite physical optimization and exact rational endpoint replay only',
            'classes':35,'starts_per_class':4,'completed_runs':len(runs),
            'optimizer_successes':sum(r['success']for r in runs),'largest_iterations':max(r['iterations']for r in runs),
            'approximate_minimum':min(r['value']for r in runs),'approximate_maximum_at_stops':max(r['value']for r in runs),
            'max_unitarity_error':max(r['unitarity_max_error']for r in runs),
            'gradient_max_abs_error':search['gradient_max_abs_error'],'exact_endpoints':len(values),
            'exact_negative':sum(v<0 for v in values),'exact_zero':sum(v==0 for v in values),'exact_positive':sum(v>0 for v in values),
            'smallest_exact_endpoint':str(min(values)),
            'universal_sign_claim':False}


if __name__=='__main__':
    if '--prepare' in sys.argv:prepare()
    report=verify()
    if '--write-report' in sys.argv:(HERE/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    else:require(report==json.loads((HERE/'report.json').read_text()),'archived report')
    print(json.dumps(report,indent=2))
