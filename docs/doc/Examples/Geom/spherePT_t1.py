# - sphere (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C
import KCore.test as test

a = D.sphere((0,0,0), 1., 20)
t = C.newPyTree(['Base',2,a])
test.testT(t, 1)
test.writeCoverage(100)
