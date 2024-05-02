# - getTangent (array) -
import Geom as D
import KCore.test as test

# Along spline
c = D.polyline([(0,0,0),(1,1,0),(2,-1,0)])
a = D.spline(c, order=3, density=10.)
res = D.getTangent(a)
test.testA([res], 1)

# Along a set of splines
res =  D.getTangent([a,a,a])
test.testA(res, 2)
