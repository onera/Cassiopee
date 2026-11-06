# - blankCells (array) -
import Converter as C
import Geom as D
import Generator as G
import Transform as T
import Connector as X
import KCore.test as test

# cas 3D masque sphere
surf = D.sphere((0,0,0), 0.5, 20)
surf = T.rotate(surf,(0.,0.,0.),(0.,1.,0.), 90.)

a = G.cart((-1.,-1.,-1.),(0.1,0.1,0.1), (20,20,20))
ca = C.array('cellN',19,19,19)
ca = C.initVars(ca, 'cellN', 1.)

blankingTypes = [-1,-2,1,2]; deltas = [0.,0.1]; dim = 3; isNot = [0,1]
c = 1
for delta in deltas:
    for masknot in isNot:
        for type in blankingTypes:
            if type ==-1 and delta > 0: c += 1
            elif type ==-2 and delta > 0: c += 1
            elif type == 2 and delta > 0: c += 1
            elif masknot == 1 and delta > 0: c+= 1
            else:
                celln = X.blankCells([a], [ca], [surf], type,
                                     delta, dim, masknot)

                test.testA(celln, c)
                ca2 = C.initVars(ca, 'cellN', 1.)
                X._blankCells([a], [ca2], [surf], type,
                              delta, dim, masknot)
                test.testA([ca2], c)
                c += 1
