# - convertPyTree2File (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

LOCAL = test.getLocal()

a = G.cart((0.,0.,0.),(0.1,0.1,0.1),(11,11,11))
a = C.initVars(a,'F',1.); a = C.initVars(a,'centers:G',2.)
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
b = G.cartHexa((0.,0.,0.),(0.1,0.1,0.1),(11,11,11)); b[0] = 'cartHexa'
b = C.initVars(b,'F',1.); b = C.initVars(b,'centers:G',2.)
t = C.newPyTree(['Base',a,b])
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
C.convertPyTree2File(t, LOCAL+'/in.cgns')

# ADF ecrit la date dans les fichiers, on ne peut donc pas comparer byte a byte
#test.testF('in.cgns',1)
t = C.convertFile2PyTree(LOCAL+'/in.cgns')
test.testT(t, 1)
