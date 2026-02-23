# - adaptZoneBCEltRange2EltList (pyTree) -
import Converter.Internal as Internal
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cartHexa((0,0,0), (1,1,1), (10,10,10))
b = G.cartHexa((0,0,0), (1,1,1), (10,10,1))
a = C.addBC2Zone(a, 'wall', 'BCWall', subzone=b)

a = Internal.adaptZoneBCEltRange2EltList(a)
C.convertPyTree2File(a, 'out.cgns')
