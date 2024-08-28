# - sphere6 (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C
import KCore.test as test

A = D.sphere6((0,0,0), 1., 20)
t = C.newPyTree(['Base',2]); t[2][1][2] += A
test.testT(t, 1)
test.writeCoverage(100)
