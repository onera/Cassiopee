# - getArgMin (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T
import KCore.test as test

# test sur une zone
a = G.cart( (0,0,0), (1.,1.,1.), (10,10,10) )
a = C.initVars(a,'F',1.); a = C.initVars(a,'centers:G',2.)
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
a = T.rotate(a, (0,0,0), (0,0,1), 5.)
a = T.rotate(a, (0,0,0), (1,0,0), 5.)
a = T.rotate(a, (0,0,0), (0,1,0), 5.)
argmin = C.getArgMin(a, 'CoordinateX')
test.testO(argmin)

# Sur un arbre
b = G.cart((1., 0.2, 0.), (0.1, 0.1, 0.1), (11, 21, 2))
b = C.addBC2Zone(b, 'wall1', 'BCWall', 'imax')
t = C.newPyTree(['Base']); t[2][1][2] += [a,b]
t = C.initVars(t,'F',1.); t = C.initVars(t,'centers:G',2.)
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = T.rotate(t, (0,0,0), (0,0,1), 5.)
t = T.rotate(t, (0,0,0), (1,0,0), 5.)
t = T.rotate(t, (0,0,0), (0,1,0), 5.)
argmin = C.getArgMin(t, 'CoordinateX')
test.testO(argmin, 2)
