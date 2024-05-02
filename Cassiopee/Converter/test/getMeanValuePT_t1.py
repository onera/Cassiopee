# - getMeanValue (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

a = G.cart( (0,0,0), (1.,1.,1.), (11,1,1) )
b = G.cart( (-1,0,0), (1,1,1), (2,2,1) ); b[0] = 'cart2'
# test sur une zone
meanval = C.getMeanValue(a, 'CoordinateX')
test.testO(meanval, 1)
# test sur une liste de zones
meanval = C.getMeanValue([a,b], 'CoordinateX')
test.testO(meanval, 2)
# test sur un arbre
t = C.newPyTree(['Base',2]); t[2][1][2] += [a,b]
t = C.addBC2Zone(t, 'wall1', 'BCWall', 'imax')
t = C.initVars(t,'F',1.); t = C.initVars(t,'centers:G',2.)
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
meanval = C.getMeanValue(t, 'CoordinateX')
test.testO(meanval, 3)
