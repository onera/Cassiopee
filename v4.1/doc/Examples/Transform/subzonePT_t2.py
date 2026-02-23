# - subzone (pyTree)-
# Maillages structures + frontieres
import Converter.PyTree as C
import Transform.PyTree as T
import Generator.PyTree as G
import KCore.test as test

# 3d structure
a = G.cart((0,0,0), (1,1,1), (10,20,10))
C._addBC2Zone(a,'wall1','BCWall','jmin')
C._addBC2Zone(a, 'match1', 'BCMatch', 'imax', a, 'imin',[1,2,3])
C._addBC2Zone(a, 'overlap1', 'BCOverlap', 'jmax')
C._initVars(a, 'centers:celln',1.)
C._initVars(a,'{centers:G}={centers:CoordinateX}')
C._initVars(a,'{Density}=3*{CoordinateX}*{CoordinateY}')
t = C.newPyTree(['Base']); t[2][1][2].append(a)
t[2][1] = C.addState(t[2][1], 'EquationDimension', 3)
t = T.subzone(t, (3,1,3), (7,8,5))
test.testT(t,1)
#
# 2D structure
#
a = G.cart((0,0,0), (1,1,1), (10,20,1))
C._addBC2Zone(a,'wall1','BCWall','jmin')
C._addBC2Zone(a, 'match1', 'BCMatch', 'imax', a, 'imin',[1,2])
C._addBC2Zone(a, 'overlap1', 'BCOverlap', 'jmax')
C._initVars(a, 'centers:celln',1.)
C._initVars(a,'{centers:G}={centers:CoordinateX}')
C._initVars(a,'{Density}=3*{CoordinateX}*{CoordinateY}')
t = C.newPyTree(['Base',2]); t[2][1][2].append(a)
t[2][1] = C.addState(t[2][1], 'EquationDimension', 2)
t = T.subzone(t, (3,1,1), (7,8,1))
test.testT(t,2)

# 3d non structure
a = G.cartTetra((0,0,0), (1,1,1), (10,20,10))
C._initVars(a, 'centers:celln',1.)
C._initVars(a,'{centers:G}={centers:CoordinateX}')
C._initVars(a,'{Density}=3*{CoordinateX}*{CoordinateY}')
t = C.newPyTree(['Base']); t[2][1][2].append(a)
t[2][1] = C.addState(t[2][1], 'EquationDimension', 3)
t = T.subzone(t, [10,19,20,220])
test.testT(t,3)
#
# 2D non structure
#
a = G.cartTetra((0,0,0), (1,1,1), (10,20,1))
C._initVars(a, 'centers:G',2.)
C._initVars(a,'{F}=3*{CoordinateX}*{CoordinateY}')
t = C.newPyTree(['Base',2]); t[2][1][2].append(a)
t[2][1] = C.addState(t[2][1], 'EquationDimension', 2)
t = T.subzone(t, [8,19,18,29])
test.testT(t,4)
