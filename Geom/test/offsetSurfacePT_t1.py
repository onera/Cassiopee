# - offsetSurface (pyTree) -
import Geom.PyTree as D
import KCore.test as test

#2D
a = D.circle((0,0,0), 1)
b = D.offsetSurface(a, offset=1., pointsPerUnitLength=10., algo=0, dim=2)
test.testT(b)
b = D.offsetSurface(a, offset=1., pointsPerUnitLength=10., algo=1, dim=2)
test.testT(b,2)

# 3D
a = D.sphere((0,0,0), 1)
b = D.offsetSurface(a, offset=1., pointsPerUnitLength=10., algo=0, dim=2)
test.testT(b,3)
b = D.offsetSurface(a, offset=1., pointsPerUnitLength=10., algo=1, dim=2)
test.testT(b,4)
