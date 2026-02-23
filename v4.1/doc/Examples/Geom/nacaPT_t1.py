# - naca (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C
import KCore.test as test

a = D.naca(12.)
t = C.newPyTree(['Base',1, a])
test.testT(t, 1)

test.writeCoverage(100)
