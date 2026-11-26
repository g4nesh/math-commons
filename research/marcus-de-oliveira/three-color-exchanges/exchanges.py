#!/usr/bin/env python3
"""Exact local-exchange kernel checks; no optimizer or symbolic-algebra dependency."""
from fractions import Fraction as F
from itertools import combinations,permutations,product
from pathlib import Path
import importlib.util,sys,json
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('exchange_coupling',HERE.parent/'three-color-coupling'/'coupling.py')
c=importlib.util.module_from_spec(spec);sys.modules[spec.name]=c;spec.loader.exec_module(c)


def colorings(pattern):
    if not pattern or any(type(r)is not int or r<0 for r in pattern):raise ValueError('Expected nonnegative integer color counts')
    n=sum(pattern)
    if n>7:raise ValueError('Finite replay is bounded to dimension7; the theorem has no such limit')
    return tuple(word for word in product([j for j,r in enumerate(pattern)if r],repeat=n)if tuple(word.count(j)for j in range(len(pattern)))==tuple(pattern))


def statistic_rows(words,color_count):
    n=len(words[0]);out=[]
    for k in range(n+1):
        for subset in combinations(range(n),k):
            counts=[tuple(sum(word[i]==j for i in subset)for j in range(color_count))for word in words]
            for count in sorted(set(counts)):out.append([int(x==count)for x in counts])
    return out


def moves(pattern,words=None):
    words=colorings(pattern)if words is None else words
    active=[j for j,r in enumerate(pattern)if r]
    if len(active)<3:return [],[]
    n=sum(pattern);index={word:i for i,word in enumerate(words)};rows=[];labels=[]
    for sites in combinations(range(n),3):
        other=[i for i in range(n)if i not in sites]
        for selected_colors in combinations(active,3):
            residual=tuple(r-int(j in selected_colors)for j,r in enumerate(pattern))
            for outside in colorings(residual):
                vector=[0]*len(words)
                for colors in permutations(selected_colors):
                    word=[None]*n
                    for i,j in zip(sites,colors):word[i]=j
                    for i,j in zip(other,outside):word[i]=j
                    sign=(-1)**sum(colors[i]>colors[j]for i in range(3)for j in range(i+1,3))
                    vector[index[tuple(word)]]=sign
                rows.append(vector);labels.append({'sites':list(sites),'colors':list(selected_colors),'outside':list(outside)})
    return rows,labels


def rank_and_rows(rows):
    if not rows:return 0,[]
    selected=c.independent_rows(rows)
    return len(selected),selected


def solve_coefficients(basis,vector):
    if any(isinstance(v,bool)or not isinstance(v,(int,F))for v in vector) or any(any(isinstance(v,bool)or not isinstance(v,(int,F))for v in row)for row in basis):
        raise TypeError('Exchange solving requires exact integer or Fraction scalars')
    if any(len(row)!=len(vector)for row in basis):raise ValueError('Wrong exchange-vector dimensions')
    if not basis:
        if any(vector):raise ValueError('Nonzero vector with empty kernel basis')
        return []
    a=[[F(row[i])for row in basis]+[F(vector[i])]for i in range(len(vector))];n=len(basis);pivot=0;columns=[]
    for j in range(n):
        chosen=next((i for i in range(pivot,len(a))if a[i][j]),None)
        if chosen is None:raise ValueError('Dependent exchange basis')
        a[pivot],a[chosen]=a[chosen],a[pivot];factor=a[pivot][j];a[pivot]=[x/factor for x in a[pivot]]
        for i in range(len(a)):
            if i!=pivot and a[i][j]:
                factor=a[i][j];a[i]=[x-factor*y for x,y in zip(a[i],a[pivot])]
        columns.append(j);pivot+=1
    if any(all(x==0 for x in row[:-1])and row[-1]for row in a):raise ValueError('Vector is outside the exchange span')
    result=[a[i][-1]for i in range(n)]
    if any(sum(coeff*row[j]for coeff,row in zip(result,basis))!=F(vector[j])for j in range(len(vector))):raise ValueError('Exchange reconstruction failed')
    return result


def check_pattern(pattern):
    words=colorings(pattern);A=statistic_rows(words,len(pattern));H,labels=moves(pattern,words)
    if any(sum(x*y for x,y in zip(a,h))for a in A for h in H):raise ValueError('A move changes a grouped statistic')
    rank_a,_=rank_and_rows(A);rank_h,selected=rank_and_rows(H)
    if rank_h!=len(words)-rank_a:raise ValueError('Local moves fail to span the kernel')
    return {'pattern':list(pattern),'colorings':len(words),'statisticRows':len(A),'statisticRank':rank_a,
            'localMoves':len(H),'exchangeRank':rank_h,'kernelDimension':len(words)-rank_a},(H,labels,selected)


def generate():
    report,(H,labels,selected)=check_pattern((2,2,2));data=c.load_data();correction=list(map(F,data['correction']))
    coefficients=solve_coefficients([H[i]for i in selected],correction)
    return {'pattern':[2,2,2],'basisMoveIndices':selected,'basisMoveLabels':[labels[i]for i in selected],
            'interiorRepairCoefficients':list(map(str,coefficients))}


def four_color_dual():
    words=colorings((1,1,1,1));H,_=moves((1,1,1,1),words)
    nums=[[-31,492,-240,300],[-300,225,500,0],[240,-180,225,500],[-492,-256,-180,225]]
    u=[[c.Q(F(v,625))for v in row]for row in nums]
    canonical=[]
    for word in words:
        z=c.Q((-1)**sum(word[i]>word[j]for i in range(4)for j in range(i+1,4)))*c.det(u).conjugate()
        for i,j in enumerate(word):z*=u[i][j]
        canonical.append(z.real)
    y=[F(sum(tuple(sorted(word[i]for i in rows))==columns for rows,columns in c.COVER[:13])-1)for word in words]
    if min(y)<0 or any(sum(a*b for a,b in zip(y,h))for h in H):raise ValueError('Invalid four-color cone dual')
    gap=sum(a*b for a,b in zip(y,canonical))
    if gap!=F(-672,15625):raise ValueError('Four-color separating gap changed')
    return {'pattern':[1,1,1,1],'nonnegativeDualEntries':len(y),'exchangeDirectionsAnnihilated':len(H),'canonicalDualValue':str(gap)}


def verify(data):
    report,(H,labels,selected)=check_pattern((2,2,2))
    if data['pattern']!=[2,2,2]or data['basisMoveIndices']!=selected or data['basisMoveLabels']!=[labels[i]for i in selected]:raise ValueError('Exchange basis changed')
    if len(data['interiorRepairCoefficients'])!=47 or any(type(x)is not str for x in data['interiorRepairCoefficients']):raise ValueError('Expected47 exact rational coefficients')
    coefficients=list(map(F,data['interiorRepairCoefficients']));source=c.load_data();correction=list(map(F,source['correction']))
    if any(sum(x*H[k][j]for x,k in zip(coefficients,selected))!=correction[j]for j in range(90)):raise ValueError('Published interior repair is not reconstructed')
    # Gauge rotations alter bases within each color eigenspace, not the weights.
    u=c.fixture();g=c.e.identity(6)
    for j,t in enumerate((F(1,2),F(1,3),F(1,4))):
        g=c.e.multiply(c.rotation(6,2*j,2*j+1,t,c.Q(0,1)if j%2 else c.O),g)
    rotated=c.e.multiply(u,g);c.check_unitary(rotated)
    if c.canon(rotated)!=c.canon(u):raise ValueError('Within-color gauge invariance failed')
    small=[]
    for pattern in ((0,0,0),(0,2,2),(1,1,1),(2,1,1),(2,2,1),(1,1,1,1)):
        case,_=check_pattern(pattern);small.append(case)
    return {'main':report,'degenerateAndSmallerPatterns':small,'interiorRepairDirections':47,
            'interiorRepairCoordinatesChecked':90,'gaugeRotationsChecked':3,'unchangedNegativeCanonicalMasses':sum(x<0 for x in c.canon(u)),
            'normalizedDualPolytopeDimension':42,'knownFourColorDual':four_color_dual(),'scope':'Complete linear repair space and cone-dual reduction. No universal positivity theorem or Conjecture F counterexample.'}


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--derive',action='store_true');parser.add_argument('--report',type=Path)
    args=parser.parse_args();path=HERE/'basis.json'
    if args.derive:path.write_text(json.dumps(generate(),indent=2)+'\n')
    result=verify(json.loads(path.read_text()))
    if args.report:args.report.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
