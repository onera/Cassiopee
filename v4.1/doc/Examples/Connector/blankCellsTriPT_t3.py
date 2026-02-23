# - blankCellsTri (array) - 'NODE IN'
import Converter.PyTree as C
import Connector.PyTree as X
import Generator.PyTree as G
import Geom.PyTree as D
import KCore.test as test
import Transform.PyTree as T

# Test 1
# Tet mask

s1 = D.sphere6((0,0,0), 1., 20, ntype='TRI')
s2 = D.sphere6((1,0,0), 1.5, N=20, ntype='TRI')
s2 = T.reorder(s2, (1,))
#C.convertPyTree2File(s2, 's2.cgns')

snear = 0.15
hWall = 0.1
raison = 1.2
smoothIter = 20
nlayer = 3
dlayer = hWall*(1.-raison**nlayer)/(1.-raison);
d = G.cart((0.,0.,0.), (dlayer/nlayer,1,1), (nlayer,1,1))
penta = G.addNormalLayers(s2, d, check=0, niter=smoothIter)
#C.convertPyTree2File(penta, 'penta.cgns')
t = C.newPyTree(['Pentas']); t[2][1][2].append(penta)
# Blanking
C._initVars(t, 'centers:cellN', 1.)
t = X.blankCellsTri(t, [[s1]], [], blankingType="cell_intersect", tol=1.e-12)
#C.convertPyTree2File(t, 'out.cgns')
test.testT(t,1)

s2 = D.sphere6((1,0,0), 1.5, N=20, ntype='QUAD')
pyras = G.quad2Pyra(s2, hratio=1.)
t = C.newPyTree(['Pyras']); t[2][1][2].append(pyras)
# Blanking
C._initVars(t, 'centers:cellN', 1.)
t = X.blankCellsTri(t, [[s1]], [], blankingType="cell_intersect", tol=1.e-12)
#C.convertPyTree2File(t, 'out2.cgns')
test.testT(t,2)
