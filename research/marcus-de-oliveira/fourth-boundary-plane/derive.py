#!/usr/bin/env python3
"""Optional exact Cramer interpolation; no symbolic-algebra or LP package."""
import argparse,json
from fractions import Fraction as F
from pathlib import Path
import plane as p


def interpolate(grid):
    xcoefs=[p.m.interpolate(row)for row in grid];out={}
    for i in range(13):
        coefficients=p.m.interpolate([row[i]if i<len(row)else F(0)for row in xcoefs])
        for j,v in enumerate(coefficients):
            if v:out[i,j]=v
    if any(i+j>12 for i,j in out):raise ValueError('Determinant degree bound failed')
    return out


def encode(poly):return [[i,j,str(v)]for(i,j),v in sorted(poly.items())]

def derive():
    grids=[[[]for _ in range(13)]for _ in range(37)]
    for t in range(13):
        for s in range(13):
            for grid,v in zip(grids,p.raw_values(F(s),F(t))):grid[t].append(v)
    polys=[interpolate(grid)for grid in grids]
    return {'coordinates':'a3=s*(1-i)+t*(1+i)','support':[list(x)for x in p.SUPPORT],
            'basisRows':p.ROWS,'orientation':1,'denominator':encode(polys[0]),
            'numerators':[encode(v)for v in polys[1:13]],'dualColumns':[encode(v)for v in polys[13:]],
            'sInterval':list(map(str,p.SI)),'tInterval':list(map(str,p.TI)),'boundaryPolynomial':encode(p.Q)}

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();data=derive();args.output.write_text(json.dumps(data,indent=2)+'\n')
    print('Reconstructed37 exact bivariate determinant polynomials.')
