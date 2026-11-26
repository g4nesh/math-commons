"""Exact physical failure of one parent shared by two row orientations.

Python 3.9+ stdlib. Includes positive laws for both orientations separately.
"""
from fractions import Fraction as F
from itertools import combinations
from math import prod
from pathlib import Path
import importlib.util
import json
import sys

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('parent_orientation', HERE.parent / 'split-row-orientations' / 'check.py')
orientation = importlib.util.module_from_spec(spec)
spec.loader.exec_module(orientation)
base = orientation.base
Q, exact, require, dot = base.Q, base.exact, base.require, base.dot


def load():
    return json.loads((HERE / 'certificate.json').read_text())


def response(u):
    """The PR102 symbolic Pauli expansion, at this different coarse unitary."""
    base.check_unitary(u)
    det = exact.determinant(u)
    s = base.structure()
    out = [[Q(0)]*4 for _ in range(33)]
    for p in base.PERMS:
        factor = base.parity(p)*prod(u[i][p[i]] for i in range(4))/det
        a,b,c,d = u[4][p[4]],u[5][p[4]],u[4][p[5]],u[5][p[5]]
        terms = ((a*d-b*c)/2,(-a*c+b*d)/2,Q(0,-1)*(a*c+b*d)/2,(a*d+b*c)/2)
        row = out[s['lookup'][base.table(p)]]
        for j,value in enumerate(terms):
            row[j] += factor*value
    return [[x.real for x in row] for row in out]


def physical_coefficients(frame, norm):
    """All 924 squared minors of diag(I4,norm^-1/2 I2)*frame."""
    s = base.structure()
    result = {key:F(0) for key in s['keys']}
    count = 0
    for size in range(7):
        for rows in combinations(range(6),size):
            for cols in combinations(range(6),size):
                value = exact.determinant([[frame[i][j] for j in cols] for i in rows])
                squared = (value.real**2+value.imag**2)/norm**sum(i>=4 for i in rows)
                key = tuple(sum(i not in rows and base.ROW[i]==r for i in range(6)) for r in range(4))
                key += tuple(sum(base.COL[j]==c for j in cols) for c in range(3))
                result[key] += squared
                count += 1
    require(count==924,'all squared minors')
    return result


def verify(data=None, minors=True):
    data = load() if data is None else data
    require(data['schema']=='mdo.common-parent-obstruction.v1','schema')
    base.verify(minors=False)
    s = base.structure()
    u = base.fixture(data)
    r = response(u)
    require(all(sum(r[j][k] for j in f)==0 for f in s['fibers'] for k in (1,2,3)), 'invariant coarse parents')
    require(all(r[f[0]][0]==r[f[1]][0] for _,f in s['doubles']), 'odd constant')
    coarse = [sum(r[j][0] for j in f) for f in s['fibers']]
    d = [[r[f[0]][k]-r[f[1]][k] for k in (1,2,3)] for _,f in s['doubles']]
    t = [[sum(row[j]*d[j][k] for j in range(12)) for k in range(3)] for row in s['b']]
    require(coarse==list(map(F,data['parent_canonical'])),'coarse certificate')
    require(t==[list(map(F,row)) for row in data['moment_response']],'moment certificate')
    # Parent positivity rows followed by the31 support inequalities A*r <= c-||T^T h||.
    a = [[-s['kernel'][v][p] for v in range(6)] for p in range(21)]
    constants, norm_squares, magnitudes = [], [], []
    for h in s['normals']:
        mag = [abs(sum(h[i]*s['b'][i][j] for i in range(6))) for j in range(12)]
        magnitudes.append(mag)
        a.append([-sum(m*s['kernel'][v][p] for m,(p,_) in zip(mag,s['doubles'])) for v in range(6)])
        constants.append(sum(m*coarse[p] for m,(p,_) in zip(mag,s['doubles'])))
        norm_squares.append(sum(sum(h[i]*t[i][k] for i in range(6))**2 for k in range(3)))
    dual = list(map(F,data['dual']))
    require(len(dual)==52 and all(x>=0 for x in dual),'nonnegative52-row dual')
    require([(i,x) for i,x in enumerate(dual) if x]==[(30,F(1,2)),(47,F(1,2))],'two-normal support')
    require(all(sum(dual[i]*a[i][j] for i in range(52))==0 for j in range(6)),'all six repairs cancel')
    lower = list(map(F,data['norm_lower_bounds']))
    require(len(lower)==31 and all(x>=0 and x*x<=v for x,v in zip(lower,norm_squares)), 'certified norm lower bounds')
    bound = dot(dual,coarse+[c-l for c,l in zip(constants,lower)])
    require(bound==F(data['dual_upper_bound'])<0,'negative exact sphere dual bound')
    require(len(data['directions'])==2,'two physical directions')
    responses, floors, parents, directions = [], [], [], []
    for record,expected in zip(data['directions'],(9,26)):
        require(record['normal_index']==expected,'normal order')
        n = orientation.rational_direction(tuple(map(F,record['direction'])))
        q = [row[0]+dot(row[1:],n) for row in r]
        rotation = orientation.direction_frame(n)
        norm = 2*(1+n[2]) if n[2]!=-1 else F(1)
        require(exact.multiply(exact.adjoint(rotation),rotation)==[[Q(norm),Q(0)],[Q(0),Q(norm)]], 'normalized physical frame')
        frame = orientation.rotate(rotation,u)
        require(orientation.canonical_frame(frame)==q,'direct grouped determinant identity')
        weights = list(map(F,record['weights']))
        require(len(weights)==33 and sum(weights)==1 and min(weights)==F(record['floor'])>0,'individual strict law')
        require(all(dot(row,[w-v for w,v in zip(weights,q)])==0 for row in s['matrix']),'all192 individual polynomial identities')
        if minors:
            physical = physical_coefficients(frame,norm)
            require(all(dot(row,weights)==physical[key] for row,key in zip(s['matrix'],s['keys'])),'independent924-minor physical law')
        h = s['normals'][expected]
        value = dot([sum(h[i]*t[i][k] for i in range(6)) for k in range(3)],n)
        require(value==F(record['response']),'rational directional action')
        responses.append(value);floors.append(min(weights));directions.append(list(map(str,n)))
        parents.append([sum(weights[j] for j in f) for f in s['fibers']])
    require(parents[0]!=parents[1],'the two exhibited parents differ')
    pair_bound = (constants[9]+constants[26]-sum(responses))/2
    require(pair_bound==F(data['two_direction_dual_upper_bound'])<0,'two directions already obstruct one parent')
    # A positive open neighborhood retains both strict laws and the incompatibility.
    radius = F(1,100000000)
    lipschitz = 1344*(sum(magnitudes[9])+sum(magnitudes[26]))
    require(lipschitz*radius < -pair_bound/2,'strict obstruction neighborhood')
    require(1344*radius < min(floors)/2,'strict individual law neighborhood')
    return {'physical_coarse_fixture':'four exact complex Householder reflections',
            'normal_indices':[9,26], 'normals':[list(s['normals'][i]) for i in (9,26)],
            'six_repair_rows':[a[21+i] for i in (9,26)],
            'sphere_dual_upper_bound':str(bound),'two_direction_dual_upper_bound':str(pair_bound),
            'common_parent_for_these_two_directions':False,'separate_positive_laws_for_both':True,
            'individual_mass_floors':list(map(str,floors)), 'coarse_entrywise_unitary_radius':str(radius),
            'obstruction_lipschitz_constant':lipschitz,'individual_lipschitz_constant':1344,
            'nearby_individual_floor':str(min(floors)/2),
            'polynomial_coefficients_checked':384,'squared_minors_checked':1848 if minors else 0,
            'zero_entries':sum(x==Q(0) for row in u for x in row),
            'nonreal_entries':sum(x.imag!=0 for row in u for x in row),
            'discovery_cases_screened':4, 'discovery_seed':data['seed']}


if __name__=='__main__':
    report = verify()
    if '--write-report' in sys.argv:
        (HERE/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    else:
        require(report==json.loads((HERE/'report.json').read_text()),'archived report')
    print(json.dumps(report,indent=2))
