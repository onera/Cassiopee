# - adaptBCFace2BCC (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G

a = G.cartHexa((0,0,0), (1,1,1), (10,10,10))
a = C.addBC2Zone(a, 'wall', 'BCWall', faceList=[1,7])
a = Internal.adaptBCFace2BCC(a)
C.convertPyTree2File(a, 'out.cgns')
