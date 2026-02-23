# - splitCurvatureAngle (pyTree)-
import Converter.PyTree as C
import Geom.PyTree as D
import Transform.PyTree as T
import KCore.test as test

# i-array
a = D.naca(12, 5001)
a = C.initVars(a,'Density',2.); a = C.initVars(a,'centers:cellN',1.)
zones = T.splitCurvatureAngle(a, 20.)
test.testT(zones)
