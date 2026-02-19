# - adaptBCFacePL2VertexPL (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G
import KCore.test as test

# - NGons -

t = G.cartNGon((2,0,0), (0.1,0.1,1), (10,10,2))
t = C.addBC2Zone(t, 'wall', 'BCWall', faceList=[1, 11, 21, 31, 41])
t = C.addBC2Zone(t, 'wall', 'BCViscousWall', faceList=[10, 20, 30, 40, 50])
C._initBCDataSet(t, '{var}=1.')

# Adapt face PL into vertex PL for all BC nodes, remove BCDatasets
a = Internal.copyTree(t)
a = Internal.adaptBCFacePL2VertexPL(a, remove=True)
test.testT(a, 1)

# Adapt face PL into vertex PL for all BC nodes with name wall* (in place), keep BCDatasets
a = Internal.copyTree(t)
bcs = Internal.getNodesFromName(a, 'wall*')
Internal._adaptBCFacePL2VertexPL(a, bcs=bcs, remove=False)
test.testT(a, 2)

# Adapt face PL into vertex PL for all BCViscousWall nodes (in place), keep BCDatasets
a = Internal.copyTree(t)
Internal._adaptBCFacePL2VertexPL(a, btype='BCViscousWall', remove=False)
test.testT(a, 3)
