# - adapt2FastP (pyTree) -
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.PyTree as C
import Converter.Internal as Internal
import Connector.PyTree as X
import Converter.GhostCells as GC
import KCore.test as test

a = G.cartNGon((0,0,0), (1,1,1), (6,6,3))
b = T.splitNParts(a, N=3)
t = C.newPyTree(['Base',b])
t = X.connectMatch(t)
t = C.fillEmptyBCWith(t, 'wall', 'BCWall', dim=3)
Internal._adaptNFace2PE(t, remove=False) 
# Test avec deux couches
t = GC.adapt2FastP(t, nlayers=2)       # creation parentElement du NGon
test.testT(t, 1)

