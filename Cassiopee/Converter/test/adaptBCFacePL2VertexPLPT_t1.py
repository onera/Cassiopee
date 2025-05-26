# - adaptBCFacePL2VertexPL (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G
import KCore.test as test

# - NGons -

# Adapt face PL into vertex PL for all BC nodes, remove BCDatasets
a = G.cartNGon((2,0,0), (0.1,0.1,1), (10,10,2))
a = C.addBC2Zone(a, 'wall', 'BCWall', faceList=[1,2])
a = C.addBC2Zone(a, 'wall', 'BCViscousWall', faceList=[3,5])
C._initBCDataSet(a, '{var}=1.')
a = Internal.adaptBCFacePL2VertexPL(a, remove=True)
test.testT(a, 1)

# Adapt face PL into vertex PL for all BC nodes with name wall* (in place), keep BCDatasets
a = G.cartNGon((2,0,0), (0.1,0.1,1), (10,10,2))
a = C.addBC2Zone(a, 'wall', 'BCWall', faceList=[1,2])
a = C.addBC2Zone(a, 'wall', 'BCViscousWall', faceList=[3,5])
C._initBCDataSet(a, '{var}=1.')
bcs = Internal.getNodesFromName(a, 'wall*')
Internal._adaptBCFacePL2VertexPL(a, bcs=bcs, remove=False)
test.testT(a, 2)

# Adapt face PL into vertex PL for all BCViscousWall nodes (in place), keep BCDatasets
a = G.cartNGon((2,0,0), (0.1,0.1,1), (10,10,2))
a = C.addBC2Zone(a, 'wall', 'BCWall', faceList=[1,2])
a = C.addBC2Zone(a, 'wall', 'BCViscousWall', faceList=[3,5])
C._initBCDataSet(a, '{var}=1.')
Internal._adaptBCFacePL2VertexPL(a, btype='BCViscousWall', remove=False)
test.testT(a, 3)
