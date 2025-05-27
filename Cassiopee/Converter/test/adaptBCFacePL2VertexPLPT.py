# - adaptBCFacePL2VertexPL (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G

# - NGons -

# Adapt face PL into vertex PL for all BC nodes with name 'wall*', keep BCDatasets
a = G.cartNGon((2,0,0), (0.1,0.1,1), (10,10,2))
a = C.addBC2Zone(a, 'wall', 'BCWall', faceList=[1,2])
a = C.addBC2Zone(a, 'wall', 'BCViscousWall', faceList=[3,5])
C._initBCDataSet(a, '{var}=1.')
bcs = Internal.getNodesFromName(a, 'wall*')
Internal._adaptBCFacePL2VertexPL(a, bcs=bcs, remove=False)
C.convertPyTree2File(a, 'out.cgns')
