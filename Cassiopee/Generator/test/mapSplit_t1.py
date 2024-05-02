# - mapSplit (array) -
import Generator as G
import Geom as D
import KCore.test as test

# polyline
a = D.polyline([(0.,0.,0.),(1.,1.,0.),(2.,0.,0.),(2.5,0.5,0.), (3.5,0.,0.)])

# distribution
Ni = 41
dist = G.cart((0,0,0),(1./(Ni-1),1,1),(Ni,1,1))
dist = G.enforceX(dist, 15.5/(Ni-1), 0.005, 2,5)
dist = G.enforceX(dist, 27.5/(Ni-1), 0.005, 2,5)
a = G.mapSplit(a,dist,1.e-01)
test.testA(a)
