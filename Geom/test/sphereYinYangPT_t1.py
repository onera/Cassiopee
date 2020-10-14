# - sphereYinYang (pyTree) -
import Geom.PyTree as D
import KCore.test as test

a = D.sphereYinYang((0,0,0), 1., 50)
test.testT(a, 1)
