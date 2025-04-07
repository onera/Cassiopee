# - blankCellsTetra (array) - 'NODE IN'
import Converter as C
import Connector as X
import Generator as G
import KCore.test as test

# Test 1
# Tet mask
mT4 = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
mT4 = C.convertArray2Tetra(mT4)
# Mesh to blank
a = G.cart((-5.,-5.,-5.), (0.5,0.5,0.5), (100,100,100))
# celln init
ca = C.array('cellN',100,100,100)
ca = C.initVars(ca, 'cellN', 1.)
# Blanking
celln = X.blankCellsTetra([a], [ca], mT4, blankingType=0, tol=1.e-12)
celln = C.addVars([[a], celln])
#C.convertArrays2File(celln, 'out0.plt')
test.testA(celln,1)
