# - getBCNodesFromType (pyTree) -
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
nodes = Internal.getBCNodesFromType(t, 'BCOutflow')
print([n[0] for n in nodes])
#>> ['outlet']

nodes = Internal.getBCNodesFromType(t, ['BCInflow', 'BCOutflow'])
print([n[0] for n in nodes])
#>> ['inlet', 'outlet']

nodes = Internal.getBCNodesFromType(t, 'BCWall*')
print([n[0] for n in nodes])
#>> ['wall1', 'wall2']

# On a BC_t that is not FamilySpecified
nodes = Internal.getBCNodesFromType(bc, 'BCInflow')
print([n[0] for n in nodes])
#>> ['inlet']
