# - cart (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,11,12))
b = G.cart((0.,0.,0.), (0.1,0.1,1.), (20,10,1))
c = G.cart((0.,0.,0.), (0.1,0.1,1.), (20,1,1))
t = C.newPyTree(['Base1',3,'Base2',2,'Base3',1])
t[2][1][2] += [a]; t[2][2][2] += [b]; t[2][3][2] += [c]
test.testT(t, 1)
test.writeCoverage(100)
