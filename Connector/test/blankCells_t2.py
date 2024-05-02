# - blankCells (array) -
# test pierce points XRAY
# cas surface 2D avec body en BAR
import Converter as C
import Connector as X
import Generator as G
import Geom as D
import KCore.test as test

surf = D.circle((0,0,0), 0.5, 0., 360.)
surf = C.convertArray2Tetra(surf)
res = [surf]

a = G.cart((-1.,-1.,0.),(0.1,0.1,1.), (20,20,1))
ca = C.array('cellN',19,19,1)
ca = C.initVars(ca, 'cellN', 1.)
blankingTypes = [-1,-2,1,2]; deltas = [0.,0.1]; dim = 2; isNot = [0,1]
c = 1
for delta in deltas:
    for masknot in isNot:
        for type in blankingTypes:
            celln = X.blankCells([a],[ca],[surf],type, delta, dim, masknot)
            test.testA(celln,c)
            c += 1
