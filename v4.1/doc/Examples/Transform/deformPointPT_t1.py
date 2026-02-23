# - deformPoint (PyTree) -
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.PyTree as C
import KCore.test as test

# zone structuree 2d
a = G.cart((0,0,0), (1,1,1), (10,10,2))
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
a = C.addBC2Zone(a, 'overlap1', 'BCOverlap', 'jmin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'imax', a, 'imin', [1,2])
C._addVars(a, 'F'); C._addVars(a, 'centers:G')
a = T.deformPoint(a, (0,0,0), (0.1,0.1,0.), 0.5, 0.4)
test.testT(a, 1)

# zone structuree 3d
a = G.cart((0,0,0), (1,1,1), (10,10,10))
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
a = C.addBC2Zone(a, 'overlap1', 'BCOverlap', 'jmin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'imax', a, 'imin', [1,2,3])
C._addVars(a, 'F'); C._addVars(a, 'centers:G')
a = T.deformPoint(a, (0,0,0), (0.1,0.1,0.), 0.5, 0.4)
test.testT(a,2)

# traitement sur l'arbre
a = G.cart((0,0,0), (1,1,1), (10,10,2))
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
a = C.addBC2Zone(a, 'overlap1', 'BCOverlap', 'jmin')
C._addVars(a, 'F'); C._initVars(a, 'centers:G',1.)
t = C.newPyTree(['Base',2,a])
t[2][1] = C.addState(t[2][1], 'EquationDimension', 2)
t = T.deformPoint(t, (0,0,0), (0.5,0.5,0.), 0.5, 0.4)
test.testT(t,3)
