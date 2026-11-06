# - getNearestPointIndex (array) -
import Generator as G
import Geom as D
import KCore.test as test

a = G.cart((0.,0.,0.), (0.1,0.1,0.2),(10,10,1))
inds = D.getNearestPointIndex(a, (0.54,0.34,0))
test.testO(inds, 1)
inds = D.getNearestPointIndex(a, [(0.54,0.34,0), (1,1,1)])
test.testO(inds, 2)
a2 = G.cart((2.,0.,0.), (0.1,0.1,0.2),(10,10,1))
inds = D.getNearestPointIndex([a,a2], [(0.54,0.34,0), (1,1,1)])
test.testO(inds, 3)
