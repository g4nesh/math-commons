import copy,json,unittest
from fractions import Fraction as F
import coupling as c

class CouplingTests(unittest.TestCase):
    def test_full_solver_free_replay_matches_report(self):
        self.assertEqual(c.verify(),json.loads((c.HERE/'verification.json').read_text()))

    def test_nonconstant_neighbor_weights_and_direct_grouped_minors(self):
        data=c.load_data();base=c.fixture();radius=F(c.verify(data)['closedEntrywiseUnitaryBallRadius'])
        u=c.e.multiply(c.rotation(6,0,5,radius/10,c.Q(0,1)),base)
        weights=c.corrected_weights(u,data)
        self.assertNotEqual(weights,list(map(F,data['weights'])))
        keys,rows=c.incidence()
        for subset in ((0,),(0,3),(0,2,5),(1,2,3,5)):
            for key,row in zip(keys,rows):
                if key[0]!=subset:continue
                target=sum(c.abs2(c.minor(u,subset,t))for t in c.combinations(range(6),len(subset))if tuple(sum(j in block for j in t)for block in c.BLOCKS)==key[1])
                self.assertEqual(sum(a*b for a,b in zip(row,weights)),target)
        self.assertGreaterEqual(1-c.cover_sum(u),F(c.verify(data)['obstructionGapLowerBoundThroughoutBall']))

    def test_original_determinant_with_noncollinear_three_values(self):
        u=c.fixture();weights=list(map(F,c.load_data()['weights']))
        a=[c.Q(i-2,(i*i)%5-2)for i in range(6)];beta=[c.Q(1,1),c.Q(-2,3),c.Q(4,-1)]
        b=[beta[j]for j in (0,0,1,1,2,2)]
        d=[[b[i]if i==j else c.Z for j in range(6)]for i in range(6)]
        orbit=c.e.multiply(c.e.multiply(u,d),c.e.adjoint(u))
        target=c.det([[orbit[i][j]+(a[i]if i==j else c.Z)for j in range(6)]for i in range(6)])
        predicted=c.Z
        for w,coloring in zip(weights,c.COLORINGS):
            term=c.Q(w)
            for i,color in enumerate(coloring):term*=a[i]+beta[color]
            predicted+=term
        self.assertEqual(target,predicted)

    def test_pendant_extension_all_seven_dimensional_minors(self):
        u=c.e.identity(7);base=c.fixture()
        for i in range(6):
            for j in range(6):u[i][j]=base[i][j]
        u=c.e.multiply(c.rotation(7,2,6,F(1,3),c.Q(0,1)),u)
        law={}
        for coloring,mass in zip(c.COLORINGS,map(F,c.load_data()['weights'])):
            old=tuple(coloring)+(0,);swapped=list(old);swapped[2],swapped[6]=swapped[6],swapped[2]
            for color,w in ((old,F(16,25)*mass),(tuple(swapped),F(9,25)*mass)):
                law[color]=law.get(color,F(0))+w
        colcolors=(0,0,1,1,2,2,0);checked=0
        for k in range(8):
            for subset in c.combinations(range(7),k):
                target={};predicted={}
                for columns in c.combinations(range(7),k):
                    counts=tuple(sum(colcolors[i]==j for i in columns)for j in range(3))
                    target[counts]=target.get(counts,F(0))+c.abs2(c.minor(u,subset,columns));checked+=1
                for coloring,mass in law.items():
                    counts=tuple(sum(coloring[i]==j for i in subset)for j in range(3))
                    predicted[counts]=predicted.get(counts,F(0))+mass
                self.assertTrue(all(predicted.get(key,F(0))==value for key,value in target.items()))
        self.assertEqual(checked,3432)

    def test_mass_change_is_rejected(self):
        data=c.load_data();data['weights'][0]=str(F(data['weights'][0])+1)
        with self.assertRaises(ValueError):c.verify(data)

    def test_nullspace_correction_tampering_is_rejected(self):
        data=c.load_data();data['correction'][0]=str(F(data['correction'][0])+1)
        with self.assertRaises(ValueError):c.verify(data)

    def test_color_omission_and_float_are_rejected(self):
        data=c.load_data();data['colorings'].pop()
        with self.assertRaises(ValueError):c.verify(data)
        data=c.load_data();data['weights'][0]=0.01
        with self.assertRaises(ValueError):c.verify(data)

    def test_nonunitary_and_outside_ball_rejected(self):
        u=c.fixture();u[0][0]+=1
        with self.assertRaises(ValueError):c.corrected_weights(u)
        u=c.e.multiply(c.rotation(6,0,1,F(1,100)),c.fixture())
        with self.assertRaises(ValueError):c.corrected_weights(u)

    def test_cover_permutations_and_strict_dense_fixture(self):
        u=c.fixture();self.assertTrue(all(z!=c.Z for row in u for z in row))
        self.assertEqual(len(c.COVER),21)
        self.assertGreater(1-c.cover_sum(u),F(4,100))
        self.assertEqual(sum(x<0 for x in c.canon(u)),21)

if __name__=='__main__':unittest.main()
