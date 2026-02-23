# - newPyTree (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G
import KCore.test as test

t = C.newPyTree(['Base1', 2, 'Base2', 3])
test.testT(t,1)

t = C.newPyTree(['Base1',1])
test.testT(t,2)

t = C.newPyTree(['Base1'])
test.testT(t,3)

t = C.newPyTree(['Base1','Base2'])
test.testT(t,4)

base1 = Internal.newCGNSBase('Base', 3)
base2 = Internal.newCGNSBase('Base2', 2)
t = C.newPyTree([base1,base2])
test.testT(t,5)

z1 = G.cart((0,0,0), (1,1,1), (10,10,1))
z2 = G.cart((10,0,0), (1,1,1), (10,10,1))
t = C.newPyTree(['Base', z1, z2])
test.testT(t,6)

t = C.newPyTree(['Base', z1, 'Base2', z2])
test.testT(t,7)

t = C.newPyTree(['Base', [z1,z2]])
test.testT(t,8)
