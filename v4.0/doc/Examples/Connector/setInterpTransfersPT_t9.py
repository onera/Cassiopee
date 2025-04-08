import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal
import Converter.GhostCells as GC
import Connector.PyTree as X
import KCore.test as test

N = 10; x0 = (N-1)
a = G.cartNGon((0,0,0),(1,1,1),(N,N,N))
b = G.cartNGon((x0,0,0),(1,1,1),(N,N,N))
t = C.newPyTree(['Base',a,b])
t = X.connectMatch(t)

Internal._adaptNFace2PE(t, remove=False)
t = GC.addGhostCellsNG(t, nlayers=2)
tc = C.node2Center(t)
C._initVars(tc,"{F}={CoordinateX}*{CoordinateY}")
Internal._rmNodesFromType(tc,'GridCoordinates_t')
Internal._rmNodesFromType(tc,'Elements_t')
X._setInterpData(t,tc,loc='centers',storage='inverse')
C._initVars(t,"{centers:F}=0.")
X._setInterpTransfers(t,tc,variables=['F'])
test.testT(t)
