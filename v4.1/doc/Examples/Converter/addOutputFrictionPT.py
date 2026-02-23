# - addOutputFriction (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal
import Converter.elsAProfile as elsAProfile

a = G.cart((0,0,0), (1,1,1), (10,10,10))
a = C.addBC2Zone(a, 'wall', 'BCWall', 'imin')
bcs = Internal.getNodesFromName(a, 'wall')
for b in bcs: elsAProfile._addOutputFriction(b)

C.convertPyTree2File(a, 'out.cgns')
