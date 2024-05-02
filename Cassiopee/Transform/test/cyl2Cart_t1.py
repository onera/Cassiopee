# - cyl2Cart (array) -
import Transform as T
import Generator as G
import KCore.test as test
a = G.cylinder((0.,0.,0.), 0.5, 1., 0., 360., 1., (360,20,10))
a = T.cart2Cyl(a, (0.,0.,0.),(0,0,1))
a = T.cyl2Cart(a, (0.,0.,0.),(0,0,1))
test.testA([a], 1)
