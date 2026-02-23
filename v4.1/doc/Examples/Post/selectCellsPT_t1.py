# - selectCells (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P
import KCore.test as test

def F(x, y, z):
    if (x + 2*y + z > 20.): return True
    else: return False

# CAS 1D
a = G.cart((0,0,0), (1,1,1), (30,1,1) )
C._addVars(a, 'Density'); C._addVars(a, 'centers:cellN')
a = P.selectCells(a, F, ['CoordinateX', 'CoordinateY', 'CoordinateZ'])
test.testT(a,1)

# CAS 2D
a = G.cart((0,0,0), (1,1,1), (11,11,1) )
C._addVars(a, 'Density'); C._addVars(a, 'centers:cellN')
a = C.addBC2Zone(a, 'wall', 'BCWall', 'imin')
a = P.selectCells(a, F, ['CoordinateX', 'CoordinateY', 'CoordinateZ'])
test.testT(a,2)

# CAS 3D
a = G.cart((0,0,0), (1,1,1), (11,11,11) )
C._addVars(a,'Density'); C._addVars(a, 'centers:cellN')
a = P.selectCells(a, F, ['CoordinateX', 'CoordinateY', 'CoordinateZ'])
test.testT(a,3)
