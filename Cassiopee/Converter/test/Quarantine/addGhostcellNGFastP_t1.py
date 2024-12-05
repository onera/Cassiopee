# - addGhostCellsNGFastP (pyTree) -
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.PyTree as C
import Converter.Internal as Internal
import Connector.PyTree as X
import Converter.GhostCells as GC
import KCore.test as test

a = G.cartNGon((0,0,0), (1,1,1), (100,100,3))
t = C.newPyTree(['Base',a])
C._fillEmptyBCWith(t, 'extrap', 'BCExtrapolate', dim=3)
t = T.splitNParts(t, 2)
t = X.connectMatch(t,dim=3)
Internal._adaptNFace2PE(t, remove=False)
t = GC.addGhostCellsNG(t, nlayers=2)
test.testT(t, 1)
