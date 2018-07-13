# - text2D (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C
import KCore.test as test

a = D.text2D("CASSIOPEE")
t = C.newPyTree(['Base', 2]); t[2][1][2].append(a) 
test.testT(t, 1)
test.writeCoverage(100)
