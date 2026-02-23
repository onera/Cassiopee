# - blankCells (array) -
import Converter as C
import Connector as X
import Generator as G
import Geom as D

surf = D.sphere((0,0,0), 0.5, 20)
surf = C.convertArray2Tetra(surf)

a = G.cart((-1.,-1.,-1.),(0.1,0.1,0.1), (20,20,20))
ca = C.array('cellN',19,19,19)
ca = C.initVars(ca, 'cellN', 1.)
celln = X.blankCells([a], [ca], [surf], blankingType=1, delta=0.)
a = C.node2Center(a)
celln = C.addVars([[a], celln])
C.convertArrays2File(celln, 'out.plt')

# in place-modifies cellN
surf = D.sphere((0,0,0), 0.5, 20)
surf = C.convertArray2Tetra(surf)

a = G.cart((-1.,-1.,-1.),(0.1,0.1,0.1), (20,20,20))
ca = C.array('cellN',19,19,19)
ca = C.initVars(ca, 'cellN', 1.)
X._blankCells([a], [ca], [surf], blankingType=1, delta=0.)
a = C.node2Center(a)
celln = C.addVars([[a], celln])
C.convertArrays2File(celln, 'out.plt')
