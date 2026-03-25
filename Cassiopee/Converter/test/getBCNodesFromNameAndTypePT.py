# - getBCNodesFromNameAndType (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G

a = G.cart((-0.1,0.9,0), (0.01,0.01,1.), (10,10,2))
C._addBC2Zone(a, 'inlet', 'BCInflow', 'kmin')
C._addBC2Zone(a, 'outlet', 'BCOutflow', 'kmax')
C._addBC2Zone(a, 'match1', 'BCMatch', 'imin', a, 'imax', [1,2,3])
C._addBC2Zone(a, 'wall1', 'FamilySpecified:CYLINDER', 'jmin')
C._addBC2Zone(a, 'wall2', 'FamilySpecified:CYLINDER', 'jmax')
t = C.newPyTree(['Base', a])
C._addFamily2Base(t[2][1], 'CYLINDER', bndType='BCWallViscous')
bc = Internal.getNodeFromType2(a, 'BC_t')

# On a PyTree
nodes = Internal.getBCNodesFromNameAndType(t, bndName='wall2', bndType='BCWallViscous')
print([n[0] for n in nodes])
#>> ['wall2']

nodes = Internal.getBCNodesFromNameAndType(t, bndName='wall*', bndType='BCWallViscous')
print([n[0] for n in nodes])
#>> ['wall1', 'wall2']
