# - uniformize (pyTree) -
import Geom.PyTree as D
import KCore.test as test

a = D.polyline([(0,0,0), (1,1,0), (2,0,0), (3,1,0), (4,0,0)])
a = D.uniformize(a, N=100)

test.testT(a, 1)
