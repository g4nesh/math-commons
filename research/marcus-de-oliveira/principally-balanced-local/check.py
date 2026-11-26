"""Exact fixed-basis positive law over Q(sqrt(2)); no numerical eigensolver or LP."""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations, permutations
from math import factorial
from pathlib import Path
import argparse
import json

HERE = Path(__file__).resolve().parent


@dataclass(frozen=True)
class Q2:
    a: F = F(0)
    b: F = F(0)

    def __post_init__(self):
        for name in ('a','b'):
            value=getattr(self,name)
            if isinstance(value,bool) or not isinstance(value,(int,F)):
                raise TypeError('Q2 coefficients must be exact integers or Fractions')
            object.__setattr__(self,name,F(value))

    @staticmethod
    def of(x):
        if isinstance(x, Q2):
            return x
        if isinstance(x, bool) or not isinstance(x, (int, F)):
            raise TypeError('Exact integer, Fraction, or Q2 required')
        return Q2(F(x))

    def __add__(self, x):
        x = Q2.of(x)
        return Q2(self.a+x.a, self.b+x.b)
    __radd__ = __add__

    def __neg__(self):
        return Q2(-self.a, -self.b)

    def __sub__(self, x):
        return self + -Q2.of(x)

    def __rsub__(self, x):
        return Q2.of(x) + -self

    def __mul__(self, x):
        x = Q2.of(x)
        return Q2(self.a*x.a+2*self.b*x.b, self.a*x.b+self.b*x.a)
    __rmul__ = __mul__

    def __truediv__(self, x):
        x = Q2.of(x)
        norm = x.a*x.a-2*x.b*x.b
        if norm == 0:
            raise ZeroDivisionError
        return self * Q2(x.a/norm, -x.b/norm)

    def __bool__(self):
        return bool(self.a or self.b)

    def sign(self):
        if not self.a:
            return (self.b > 0)-(self.b < 0)
        if not self.b:
            return (self.a > 0)-(self.a < 0)
        if self.a > 0 and self.b > 0:
            return 1
        if self.a < 0 and self.b < 0:
            return -1
        gap = self.a*self.a-2*self.b*self.b
        assert gap != 0
        return ((gap > 0)-(gap < 0)) * (1 if self.a > 0 else -1)

    def absolute(self):
        return self if self.sign() >= 0 else -self

    def encode(self):
        return [str(self.a), str(self.b)]


ZERO, ONE = Q2(), Q2(F(1))


def product(values):
    out = ONE
    for x in values:
        out *= x
    return out


def subsets(n):
    return [s for k in range(n+1) for s in combinations(range(n), k)]


def determinant(matrix):
    n = len(matrix)
    return sum(((-1)**sum(p[i] > p[j] for i in range(n) for j in range(i+1,n))
                * product(matrix[i][p[i]] for i in range(n))
                for p in permutations(range(n))), ZERO)


def projection_rows(n):
    return [s for k in range(n+1) for s in list(combinations(range(n), k))[1:]]


def vertex(spectrum, p, s):
    return product(spectrum[p[i]] for i in s)


def independent_columns(columns, dimension):
    echelon, chosen = [], []
    for j, col in enumerate(columns):
        v = list(col)
        for pivot, row in echelon:
            factor = v[pivot]
            v = [x-factor*y for x, y in zip(v,row)]
        pivot = next((i for i,x in enumerate(v) if x), None)
        if pivot is None:
            continue
        scale = v[pivot]
        echelon.append((pivot, [x/scale for x in v]))
        chosen.append(j)
        if len(chosen) == dimension:
            return chosen
    if len(chosen) != dimension:
        raise ValueError('Insufficient affine rank')
    return chosen


def solve(matrix, rhs):
    n = len(rhs)
    rows = [list(row)+[b] for row,b in zip(matrix,rhs)]
    for j in range(n):
        pivot = next((i for i in range(j,n) if rows[i][j]), None)
        if pivot is None:
            raise ValueError('Singular correction matrix')
        rows[j],rows[pivot] = rows[pivot],rows[j]
        scale = rows[j][j]
        rows[j] = [x/scale for x in rows[j]]
        for i in range(n):
            if i != j:
                factor = rows[i][j]
                rows[i] = [x-factor*y for x,y in zip(rows[i],rows[j])]
    return [row[-1] for row in rows]


def construct(spectrum, minors):
    spectrum = tuple(Q2.of(x) for x in spectrum)
    n = len(spectrum)
    if not n or len(set(spectrum)) != n:
        raise ValueError('Nonempty distinct real spectrum required')
    all_s = subsets(n)
    if set(minors) != set(all_s):
        raise ValueError('Every subset minor is required')
    perms = list(permutations(range(n)))
    values = [[vertex(spectrum,p,s) for s in all_s] for p in perms]
    center = {s:sum((v[i] for v in values),ZERO)/len(perms) for i,s in enumerate(all_s)}
    for k in range(n+1):
        if sum((minors[s]-center[s] for s in all_s if len(s)==k),ZERO):
            raise ValueError('Degree-sum compatibility fails')
    projected = projection_rows(n)
    indices = [all_s.index(s) for s in projected]
    cols = [[v[i]-values[0][i] for i in indices] for v in values[1:]]
    chosen = [j+1 for j in independent_columns(cols,len(projected))]
    matrix = [[values[j][i]-values[0][i] for j in chosen] for i in indices]
    correction = solve(matrix,[minors[s]-center[s] for s in projected])
    l1 = sum((x.absolute() for x in correction),ZERO)
    if (ONE/len(perms)-l1).sign() <= 0:
        raise ValueError('Sufficient strict-positivity criterion fails')
    weights = [ONE/len(perms) for _ in perms]
    for j,c in zip(chosen,correction):
        weights[j] += c
        weights[0] -= c
    assert sum(weights,ZERO) == ONE and all(w.sign() > 0 for w in weights)
    assert all(sum((w*v[i] for w,v in zip(weights,values)),ZERO)==minors[s]
               for i,s in enumerate(all_s))
    return perms,chosen,weights,l1


def transpose(a):
    return list(map(list,zip(*a)))


def multiply(a,b):
    return [[sum((x*y for x,y in zip(row,col)),F(0)) for col in zip(*b)] for row in a]


def skew_sign(n):
    return [[F((j>i)-(j<i)) for j in range(n)] for i in range(n)]


def hermitian_minors(skew):
    out = {}
    for s in subsets(len(skew)):
        det = determinant([[Q2(skew[i][j]) for j in s] for i in s])
        if len(s)%2:
            assert det == ZERO
            out[s] = ZERO
        else:
            out[s] = (-1)**(len(s)//2)*det
    return out


def report():
    # Genuine exact off-center matrix, with the known center spectrum unchanged.
    n = 4
    c,s = F(999999,1000001),F(2000,1000001)
    rotation = [[F(i==j) for j in range(n)] for i in range(n)]
    rotation[0][0]=rotation[1][1]=c
    rotation[0][1]=-s
    rotation[1][0]=s
    assert multiply(rotation,transpose(rotation)) == [[F(i==j) for j in range(n)] for i in range(n)]
    center = skew_sign(n)
    moved = multiply(multiply(rotation,center),transpose(rotation))
    assert transpose(moved)==[[-x for x in row] for row in moved]
    spectrum = [Q2(F(1),F(1)),Q2(F(-1),F(1)),Q2(F(1),F(-1)),Q2(F(-1),F(-1))]
    base_minors, minors = hermitian_minors(center),hermitian_minors(moved)
    assert len({minors[x] for x in combinations(range(4),2)})>1
    assert all(len({base_minors[x] for x in combinations(range(4),k)})==1 for k in range(5))
    for k in range(5):
        elementary = sum((product(spectrum[i] for i in a) for a in combinations(range(4),k)),ZERO)
        assert elementary == sum((base_minors[a] for a in combinations(range(4),k)),ZERO)
    perms,basis,weights,l1 = construct(spectrum,minors)
    # The general criterion, not only mere nonnegativity, has comfortable exact slack.
    assert (Q2(F(1,48))-l1).sign()>0
    assert all((w-F(1,48)).sign()>0 for w in weights)
    ranks = {}
    for size in range(1,6):
        spec = list(range(size))
        ss = subsets(size)
        center_data = {a:sum((vertex(spec,p,a) for p in permutations(range(size))),ZERO)/factorial(size) for a in ss}
        _,selected,_,_ = construct(spec,center_data)
        ranks[str(size)] = len(selected)
        assert len(selected)==2**size-size-1
    return {
        'schemaVersion':1,
        'field':'Q(sqrt(2)); [r,s] denotes r+s*sqrt(2), exact rational strings',
        'permutationConvention':'coordinate position i maps to eigenvalue label p[i], zero based',
        'center':'H0=i*K4, K4 has +1 above and -1 below diagonal',
        'offCenter':'H=R*H0*R^T; R acts on coordinates 0,1',
        'rotationCos':str(c),'rotationSin':str(s),
        'skewPartK':[[str(x) for x in row] for row in moved],
        'spectrum':[x.encode() for x in spectrum],
        'basisPermutationIndices':basis,
        'correctionL1':l1.encode(),
        'allWeightsStrictlyGreaterThan':'1/48',
        'principalMinors':[{ 'subset':list(a),'value':minors[a].encode()} for a in subsets(n)],
        'weights':[{'permutation':list(p),'weight':w.encode()} for p,w in zip(perms,weights)],
        'exactAffineRanksAtDistinctIntegerSpectra':ranks,
        'scope':'Exact off-center witness plus finite rank regressions; the all-dimensions neighborhood proof is in README.md',
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    result=report()
    path=HERE/'witness.json'
    if args.write:
        path.write_text(json.dumps(result,indent=2)+'\n')
    elif json.loads(path.read_text())!=result:
        raise SystemExit('Archived witness differs from exact replay')
    print('PASS: all 16 subset identities, 24 weights > 1/48, off-center matrix and ranks n=1..5.')
    print('Correction l1:',result['correctionL1'])


if __name__=='__main__':
    main()
