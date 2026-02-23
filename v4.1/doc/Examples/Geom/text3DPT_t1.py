# - text3D (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C
import KCore.test as test

a = D.text3D("CASSIOPEE", font='text1')
t = C.newPyTree(['Base', a])
test.testT(t, 1)
test.writeCoverage(100)
