# - distance (array) -
import Geom as D
import KCore.test as test

a = D.circle((0,0,0), 1., N=30)
b = D.circle((0,0,0), 1.1, N=30)
d = D.distance(a, b)
test.testO(d, 1)

a = D.sphere((0,0,0), 1., N=30)
b = D.sphere((0,0,0), 1.1, N=30)
d = D.distance(a, b)
test.testO(d, 2)
