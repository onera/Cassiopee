# - selectCells (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P
import KCore.test as test

def F(x, y, z):
    if (x + 2*y + z > 20.): return True
    else: return False

# CAS 1D
a = G.cart((0,9,0), (1,1,1), (11,1,1) )
C._addVars(a,'Density'); C._addVars(a,'centers:cellN')
t = C.newPyTree(['Base',1,a])
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = P.selectCells(t, F, ['CoordinateX', 'CoordinateY', 'CoordinateZ'])
test.testT(t,1)

# CAS 2D
a = G.cart((0,0,0), (1,1,1), (11,11,1) )
C._addVars(a,'Density'); C._addVars(a,'centers:cellN')
a = C.addBC2Zone(a,'ov','BCOverlap','imin')
t = C.newPyTree(['Base',2,a])
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = P.selectCells(t, F, ['CoordinateX', 'CoordinateY', 'CoordinateZ'])
test.testT(t,2)

# CAS 3D
a = G.cart((0,0,0), (1,1,1), (11,11,11) )
C._addVars(a,'Density'); C._addVars(a,'centers:cellN')
a = C.addBC2Zone(a,'ov','BCOverlap','imin')
t = C.newPyTree(['Base',3,a])
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = C.fillEmptyBCWith(t,'wall','BCWall')
t = P.selectCells(t, F, ['CoordinateX', 'CoordinateY', 'CoordinateZ'])
test.testT(t,3)
