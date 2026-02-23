# - point (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C
import KCore.test as test

a = D.point((0,0,0))
b = D.point((1,1,1))
t = C.newPyTree(['Base',1,a,b])
test.testT(t, 1)

test.writeCoverage(100)
