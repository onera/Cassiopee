# - perturbate (pyTree) -
import Transform.PyTree as T
import Generator.PyTree as G
import Converter.PyTree as C
import KCore.test as test

# Traitements sur une zone
# Structure 3D
a = G.cart( (0,0,0), (1,1,1), (10,10,10))
C._addVars(a,'Density'); C._initVars(a,'centers:cellN',1.)
a = C.addBC2Zone(a, 'wall','BCWall','imin')
a = C.addBC2Zone(a, 'overlap','BCOverlap','jmin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'imax', a, 'imin', [1,2,3])
a = T.perturbate(a, 0.1)
test.testT(a,1)

# Structure 2D
a = G.cart( (0,0,0), (1,1,1), (10,10,1))
C._addVars(a,'Density'); C._initVars(a,'centers:cellN',1.)
a = C.addBC2Zone(a, 'wall','BCWall','imin')
a = C.addBC2Zone(a, 'overlap','BCOverlap','jmin')
a = T.perturbate(a, 0.1)
test.testT(a,2)

# Traitements sur un arbre
# Structure 3D
a = G.cart( (0,0,0), (1,1,1), (10,10,10))
C._addVars(a,'Density'); C._initVars(a,'centers:cellN',1.)
a = C.addBC2Zone(a, 'wall','BCWall','imin')
a = C.addBC2Zone(a, 'overlap','BCOverlap','jmin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'imax', a, 'imin', [1,2,3])
t = C.newPyTree(['Base',a])
t[2][1] = C.addState(t[2][1], 'EquationDimension', 3)
t = T.perturbate(t, 0.1)
test.testT(t,3)

# Structure 2D
a = G.cart( (0,0,0), (1,1,1), (10,10,1))
C._addVars(a,'Density'); C._initVars(a,'centers:cellN',1.)
a = C.addBC2Zone(a, 'wall','BCWall','imin')
a = C.addBC2Zone(a, 'overlap','BCOverlap','jmin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'imax', a, 'imin', [1,2])
t = C.newPyTree(['Base',2,a])
t[2][1] = C.addState(t[2][1], 'EquationDimension', 2)
t = T.perturbate(t, 0.1)
test.testT(t,4)
