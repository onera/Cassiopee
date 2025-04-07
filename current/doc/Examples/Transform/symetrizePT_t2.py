# - symetrize (pyTree) -
# traitement sur des zones
import Transform.PyTree as T
import Generator.PyTree as G
import Converter.PyTree as C
import KCore.test as test

# Structure 3D
a = G.cart((0,0,0), (1,1,1), (10,10,3))
C._addVars(a,'Density'); C._addVars(a,'centers:cellN')
a = C.addBC2Zone(a, 'wall','BCWall','imin')
a = C.addBC2Zone(a, 'overlap','BCOverlap','jmin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'imax', a, 'imin', [1,2,3])
t = C.newPyTree(['Base',a])
t[2][1] = C.addState(t[2][1], 'EquationDimension', 3)
t = T.symetrize(t, (0.,0.,0.), (1,0,0), (0,0,1))
test.testT(t,1)

# Structure 2D
a = G.cart((0,0,0), (1,1,1), (10,10,1))
C._addVars(a,'Density'); C._addVars(a,'centers:cellN')
a = C.addBC2Zone(a, 'wall','BCWall','imin')
a = C.addBC2Zone(a, 'overlap','BCOverlap','jmin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'imax', a, 'imin', [1,2])
t = C.newPyTree(['Base',2,a])
t[2][1] = C.addState(t[2][1], 'EquationDimension', 2)
t = T.symetrize(t, (0.,0.,0.), (1,0,0), (0,0,1))
test.testT(t,2)

# TETRA
a = G.cartTetra( (0,0,0), (1,1,1), (10,10,3))
C._addVars(a,'Density'); C._addVars(a,'centers:cellN')
t = C.newPyTree(['Base',a])
t[2][1] = C.addState(t[2][1], 'EquationDimension', 3)
t = T.symetrize(t, (0.,0.,0.), (1,0,0), (0,0,1))
test.testT(t,3)

# HEXA
a = G.cartHexa((0,0,0), (1,1,1), (10,10,3))
C._addVars(a,'Density'); C._addVars(a,'centers:cellN')
t = C.newPyTree(['Base',a])
t[2][1] = C.addState(t[2][1], 'EquationDimension', 3)
t = T.symetrize(t, (0.,0.,0.), (1,0,0), (0,0,1))
test.testT(t,4)

# TRI
a = G.cartTetra((0,0,0), (1,1,1), (10,10,1))
C._addVars(a,'Density'); C._addVars(a,'centers:cellN')
t = C.newPyTree(['Base',2,a])
t[2][1] = C.addState(t[2][1], 'EquationDimension', 2)
t = T.symetrize(t, (0.,0.,0.), (1,0,0), (0,0,1))
test.testT(t,5)

# QUAD
a = G.cartHexa((0,0,0), (1,1,1), (10,10,1))
C._addVars(a,'Density'); C._addVars(a,'centers:cellN')
t = C.newPyTree(['Base',2,a])
t[2][1] = C.addState(t[2][1], 'EquationDimension', 2)
t = T.symetrize(t, (0.,0.,0.), (1,0,0), (0,0,1))
test.testT(t,6)
