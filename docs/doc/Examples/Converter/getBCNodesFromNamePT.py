# - getBCNodesFromName (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G

a = G.cart((-0.1,0.9,0), (0.01,0.01,1.), (10,10,2))
C._addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
C._addBC2Zone(a, 'match1', 'BCMatch', 'imin', a, 'imax', [1,2,3])
C._addBC2Zone(a, 'wall2', 'BCWall', 'jmax')
C._addBC2Zone(a, 'inlet', 'BCInflow', 'kmin')
C._addBC2Zone(a, 'outlet', 'BCOutflow', 'kmax')
t = C.newPyTree(['Base', a])
bc = Internal.getNodeFromType2(a, 'BC_t')

# On a PyTree
nodes = Internal.getBCNodesFromName(t, 'wall1')
print([n[0] for n in nodes])
#>> ['wall1']

nodes = Internal.getBCNodesFromName(t, ['wall1', 'outlet'])
print([n[0] for n in nodes])
#>> ['wall1', 'outlet']

nodes = Internal.getBCNodesFromName(t, 'wall*')
print([n[0] for n in nodes])
#>> ['wall1', 'wall2']

# On a BC_t
nodes = Internal.getBCNodesFromName(bc, 'wall1')
print([n[0] for n in nodes])
#>> ['wall1']
