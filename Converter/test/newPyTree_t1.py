# - newPyTree (pyTree) -
import Converter.PyTree as C
import KCore.test as test

t = C.newPyTree(['Base1', 2, 'Base2', 3])
test.testT(t,1)

t = C.newPyTree(['Base1',1])
test.testT(t,2)

t = C.newPyTree(['Base1'])
test.testT(t,3)
