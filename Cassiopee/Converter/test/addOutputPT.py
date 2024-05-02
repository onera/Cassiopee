# - addOutput(pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal
import Converter.elsAProfile as elsAProfile

a = G.cart((0,0,0), (1,1,1), (10,10,10))
C._addBC2Zone(a, 'wall', 'BCWall', 'imin')
bcs = Internal.getNodesFromName(a, 'wall')
dict = {}
dict['var'] = 'convflux_rou convflux_rov convflux_row'
dict['loc'] = 1
dict['fluxcoef'] = 1.
dict['period'] = 10
for b in bcs: elsAProfile._addOutput(b,dict)
C.convertPyTree2File(a,'out.cgns')
