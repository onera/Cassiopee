# - cyl2Cart (pyTree) -
import Transform.PyTree as T
import Generator.PyTree as G
import KCore.test as test
a = G.cylinder((0.,0.,0.), 0.5, 1., 0., 360., 1., (360,20,10))
T._cart2Cyl(a, (0.,0.,0.),(0,0,1))
T._cyl2Cart(a, (0.,0.,0.),(0,0,1))
test.testT(a, 1)
