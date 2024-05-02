# - splitCurvatureRadius (pyTree) -
import Converter.PyTree as C
import Geom.PyTree as D
import Transform.PyTree as T
import KCore.test as test

a = D.naca(12.5000)
C._addVars(a, 'Density'); C._initVars(a,'centers:cellN',1)
zones = T.splitCurvatureRadius(a, 10.)
t = C.newPyTree(['Base',1]); t[2][1][2] += zones
test.testT(t, 1)
