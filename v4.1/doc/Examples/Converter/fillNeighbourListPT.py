# - fillNeighbourList (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.elsAProfile as elsAProfile

# - Structured grids -
a = G.cylinder((0,0,0), 1., 1.5, 0., 360., 1., (80,30,2))
b = G.cart((-0.1,0.9,0), (0.01,0.01,1.), (20,20,2))
# Physical BC (here BCWall)
C._addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
# Overlap BC (with automatic donor zones)
C._addBC2Zone(a, 'overlap', 'BCOverlap', 'imin')
C._addBC2Zone(a, 'overlap1', 'BCOverlap', [1,80,30,30,1,2])
# Overlap BC (with given donor zones and doubly defined)
C._addBC2Zone(a, 'overlap2', 'BCOverlap', 'jmin', zoneDonor=[b],
              rangeDonor='doubly_defined')
C._addBC2Zone(b,"OverlapDD","BCOverlap",'imin',zoneDonor=[a],rangeDonor='doubly_defined')
t = C.newPyTree(['BaseA','BaseB'])
t[2][1][2]=[a]; t[2][2][2]=[b]
elsAProfile._overlapGC2BC(t)
elsAProfile._rmGCOverlap(t)
elsAProfile._fillNeighbourList(t)
C.convertPyTree2File(t, "out.cgns")
