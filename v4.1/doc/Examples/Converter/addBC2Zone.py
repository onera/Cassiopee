# - addBC2Zone (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

# - Structured grids -
a = G.cylinder((0,0,0), 1., 1.5, 0., 360., 1., (80,30,2))
b = G.cart((-0.1,0.9,0), (0.01,0.01,1.), (20,20,2))

# Physical BC (here BCWall)
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
# Matching BC
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'imin', a, 'imax', [1,2,3])
# Matching BC with donor zone name
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'imin', a[0], [80,80,1,30,1,2],
                 [1,2,3])
# Overlap BC (with automatic donor zones)
a = C.addBC2Zone(a, 'overlap1', 'BCOverlap', [1,80,30,30,1,2])
# Overlap BC (with given donor zones and doubly defined)
a = C.addBC2Zone(a, 'overlap2', 'BCOverlap', 'jmin', zoneDonor=[b],
                 rangeDonor='doubly_defined')
# BC defined by a family name
b = C.addBC2Zone(b, 'wall', 'FamilySpecified:myBCWall', 'imin')
# Periodic matching BC
b = C.addBC2Zone(b, 'match', 'BCMatch', 'jmin', b, 'jmax', [1,2,3],
                 translation=(0,2,0))

t = C.newPyTree(['Base',a,b])
C.convertPyTree2File(t, 'out.cgns')

# - Unstructured grids -
a = G.cartTetra((0,0,0), (1,1,1), (10,10,10))
bc = G.cartTetra((0,0,0), (1,1,1), (10,10,1))
a = C.addBC2Zone(a, 'wall1', 'BCWall', subzone=bc)
C.convertPyTree2File(a, 'out.cgns')
