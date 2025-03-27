# - blankCellsTri (array) - 'NODE IN'
import Converter as C
import Connector as X
import Generator as G
import Post as P
import Transform as T
import KCore.test as test

# Tet mask
m = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
m = P.exteriorFaces(m)
m = T.reorder(m, (-1,))
m = C.convertArray2Tetra(m)
# Mesh to blank
a = G.cart((-5.,-5.,-5.), (0.5,0.5,0.5), (100,100,100))
# celln init
ca = C.array('cellN',100,100,100)
ca = C.initVars(ca, 'cellN', 1.)
# Blanking
celln = X.blankCellsTri([a], [ca], m, blankingType=0, tol=1.e-12, cellnval=4, overwrite=1)
celln = C.addVars([[a], celln])
#C.convertArrays2File(celln, 'out0.plt')
test.testA(celln,1)
