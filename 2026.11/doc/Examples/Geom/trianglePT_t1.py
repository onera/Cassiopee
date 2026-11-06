# - triangle (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C
import KCore.test as test

a = D.triangle((0,0,0), (0.1,0.,0.1), (0.05, 0.08, 0.1))
t = C.newPyTree(['Base',2,a])
test.testT(t, 1)
test.writeCoverage(100)
