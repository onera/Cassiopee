# -fillNeigbourList (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.elsAProfile as CE
import KCore.test as test

# - Structured grids -
a = G.cylinder((0,0,0), 1., 1.5, 0., 360., 1., (80,30,2))
b = G.cart((-0.1,0.9,0), (0.01,0.01,1.), (20,20,2))

# Physical BC (here BCWall)
C._addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
# Overlap BC (with automatic donor zones)
C._addBC2Zone(a, 'overlap1', 'BCOverlap', [1,80,30,30,1,2])
# Overlap BC (with given donor zones and doubly defined)
C._addBC2Zone(a, 'overlap2', 'BCOverlap', 'jmin', zoneDonor=[b],
              rangeDonor='doubly_defined')
t = C.newPyTree(['BaseA',a,'BaseB',b])
CE._overlapGC2BC(t)
tp = CE.fillNeighbourList(t)
test.testT(tp,1)
# in place
CE._fillNeighbourList(t)
test.testT(t,1)


# - Structured grids -
a = G.cylinder((0,0,0), 1., 1.5, 0., 360., 1., (80,30,2))
b = G.cart((-0.1,0.9,0), (0.01,0.01,1.), (20,20,2))

# Physical BC (here BCWall)
C._addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
# Overlap BC (with automatic donor zones)
C._addBC2Zone(a, 'overlap1', 'BCOverlap', [1,80,30,30,1,2])
# Overlap BC (with given donor zones and doubly defined)
C._addBC2Zone(a, 'overlap2', 'BCOverlap', 'jmin', zoneDonor=['FamilySpecified:FENTE'],rangeDonor='doubly_defined')
C._tagWithFamily(b,'FENTE')
t = C.newPyTree(['Cylinder',a,'Slot',b])
C._addFamily2Base(t[2][2], 'FENTE')
CE._overlapGC2BC(t)
tp = CE.fillNeighbourList(t)
test.testT(tp,2)
# in place
CE._fillNeighbourList(t)
test.testT(t,2)
