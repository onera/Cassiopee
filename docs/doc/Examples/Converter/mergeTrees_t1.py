# - mergeTrees (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

a1 = G.cart((0.,0.,0.),(0.1,0.1,0.1),(11,11,11))
a1 = C.addBC2Zone(a1,'wall','BCWall','imin')
a1 = C.addBC2Zone(a1,'overlap','BCOverlap','imax')
a1 = C.initVars(a1,'F',1.); a1 = C.initVars(a1,'centers:G',2.)

a2 = G.cart((0.,0.,0.),(0.1,0.1,1.),(11,11,1)); a2[0] = 'cart2'
a2 = C.addBC2Zone(a2,'wall','BCWall','imin')
a2 = C.addBC2Zone(a2,'overlap','BCOverlap','imax')
a2 = C.initVars(a2,'F',1.); a2 = C.initVars(a2,'centers:G',2.)

t1 = C.newPyTree(['Base1', 2, 'Base2', 3])
t1[2][1][2].append(a2); t1[2][2][2].append(a1)
t2 = C.newPyTree(['Base2']); t2[2][1][2].append(a2)
t2[2][1] = C.addState(t2[2][1], 'Mach', 0.6)
t = C.mergeTrees(t1, t2)
test.testT(t)
