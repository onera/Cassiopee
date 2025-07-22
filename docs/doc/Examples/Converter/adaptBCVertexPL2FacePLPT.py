# - adaptBCVertexPL2FacePL (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G

# - NGons -

# Adapt vertex PL into face PL for all BC nodes with name 'wall*', keep BCDatasets
a = G.cartNGon((2,0,0), (0.1,0.1,1), (10,10,2))
a = C.addBC2Zone(a, 'wall', 'BCWall', faceList=[1, 11, 21, 31, 41])
a = C.addBC2Zone(a, 'wall', 'BCViscousWall', faceList=[10, 20, 30, 40, 50])
C._initBCDataSet(a, '{var}=1.')
bcs = Internal.getNodesFromName(a, 'wall*')
Internal._adaptBCFacePL2VertexPL(a, bcs=bcs, remove=False)  # create a vertex PL
Internal._adaptBCVertexPL2FacePL(a, bcs=bcs, remove=False)  # PL: vertex to face
C.convertPyTree2File(a, 'out.cgns')
