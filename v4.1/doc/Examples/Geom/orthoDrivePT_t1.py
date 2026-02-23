# - orthoDrive (pyTree) -
import Geom.PyTree as D
import KCore.test as test

a = D.circle((0,0,0),1.)
c = D.polyline([(0.,1.,0.), (0.,1.,1.), (2.,1.,2.)])
d = D.spline(c, 3, N=100)
o = D.orthoDrive(a, d, mode=0)
test.testT(o, 1)
