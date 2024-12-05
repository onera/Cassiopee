# - getTangent (PyTree) -
import Geom.PyTree as D
import Converter.PyTree as C
import KCore.test as test

# Along spline
c = D.polyline([(0,0,0),(1,1,0),(2,-1,0)])
a = D.spline(c, order=3, density=10.)
res =  D.getTangent(a)
test.testT(res, 1)

# Along a tree made of splines zones
t = C.newPyTree(['Base',a])
res = D.getTangent(t)
test.testT(res, 2)
