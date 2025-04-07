# - subzone (pyTree) -
# Maillages structures interieurs
import Converter.PyTree as C
import Transform.PyTree as T
import Generator.PyTree as G
import KCore.test as test

def dens(x,y): return 3*x*y

#---------------
# sur une zone
#---------------
# structure 3D + CL + champs sur une zone
a = G.cart((0,0,0), (1,1,1), (10,20,10))
C._initVars(a, 'Density', dens, ['CoordinateX','CoordinateY'])
C._initVars(a,'centers:cellN',1)
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'imax', a, 'imin',[1,2,3])
a = C.addBC2Zone(a, 'overlap1', 'BCOverlap', 'jmax')
a = T.subzone(a, (3,3,3), (7,8,5))
t = C.newPyTree(['Base',a])
test.testT(t, 1)

# structure 2D + CL
a = G.cart((0,0,0), (1,1,1), (10,20,1))
C._initVars(a,'Density',dens,['CoordinateX','CoordinateY'])
C._initVars(a,'centers:cellN',1)
a = C.addBC2Zone(a, 'wall1','BCWall', 'jmin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'imax', a, 'imin',[1,2])
a = C.addBC2Zone(a, 'overlap1', 'BCOverlap', 'jmax')
a = T.subzone(a, (3,3,1), (7,8,1))
t = C.newPyTree(['Base',2,a])
test.testT(t, 2)

#---------------
# sur un arbre
#---------------
# structure 3D + CL + champs sur une zone
a = G.cart((0,0,0), (1,1,1), (10,20,10))
C._initVars(a,'Density',dens,['CoordinateX','CoordinateY'])
C._initVars(a,'centers:cellN',1)
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'imax', a, 'imin',[1,2,3])
a = C.addBC2Zone(a, 'overlap1', 'BCOverlap', 'jmax')
t = C.newPyTree(['Base',3,a])
t[2][1] = C.addState(t[2][1], 'EquationDimension', 3)
t = T.subzone(t, (3,3,3), (7,8,5))
test.testT(t, 3)

# structure 2D + CL
a = G.cart((0,0,0), (1,1,1), (10,20,1))
C._initVars(a,'Density',dens,['CoordinateX','CoordinateY'])
C._initVars(a,'centers:cellN',1)
a = C.addBC2Zone(a, 'wall1','BCWall', 'jmin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'imax', a, 'imin',[1,2])
a = C.addBC2Zone(a, 'overlap1', 'BCOverlap', 'jmax')
t = C.newPyTree(['Base',2,a])
t[2][1] = C.addState(t[2][1], 'EquationDimension', 3)
t = T.subzone(t, (3,3,1), (7,8,1))
test.testT(t, 4)

# Indices negatifs
a = G.cart((0,0,0), (1,1,1), (10,20,10))
C._initVars(a, 'Density', dens, ['CoordinateX','CoordinateY'])
C._initVars(a,'centers:cellN',1)
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'imax', a, 'imin',[1,2,3])
a = C.addBC2Zone(a, 'overlap1', 'BCOverlap', 'jmax')
a = T.subzone(a, (3,3,3), (-2,-2,-2))
t = C.newPyTree(['Base',a])
test.testT(t, 5)
