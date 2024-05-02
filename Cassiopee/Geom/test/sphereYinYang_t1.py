# - sphereYinYang (array) -
import Geom as D
import KCore.test as test

a = D.sphereYinYang((0,0,0), 1., 50)
test.testA(a, 1)
