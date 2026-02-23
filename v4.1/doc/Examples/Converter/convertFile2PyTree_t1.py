# - convertFile2PyTree (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

LOCAL = test.getLocal()

a = G.cart((0.,0.,0.),(0.1,0.1,0.1),(11,11,11))
a = C.initVars(a,'F',1.); a = C.initVars(a,'centers:G',2.)
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
b = G.cartHexa((0.,0.,0.),(0.1,0.1,0.1),(11,11,11)); b[0] = 'cartHexa'
t = C.newPyTree(['Base']); t[2][1][2] = t[2][1][2] + [a,b]
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
C.convertPyTree2File(t, LOCAL+'/in.cgns')
t1 = C.convertFile2PyTree(LOCAL+'/in.cgns')
test.testT(t1)
