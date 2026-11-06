# - fillMissingVariables (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (10,10,11)); a[0] = 'cart1'
a = C.initVars(a,'F',1.); a = C.initVars(a,'centers:G',2.)
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
b = G.cart((1,0,0), (2,1,1), (10,10,11)); b[0] = 'cart2'
t = C.newPyTree(['Base']); t[2][1][2] += [a, b]
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = C.fillMissingVariables(t)
test.testT(t,1)
