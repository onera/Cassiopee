# - text2D (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C
import KCore.test as test

a = D.text2D("CASSIOPEE", font='text1')
t = C.newPyTree(['Base', 2, a])
test.testT(t, 1)
test.writeCoverage(100)
