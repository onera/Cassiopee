# - distrib1 (array) -
import Geom as D
import KCore.test as test

a = D.line((0,0,0), (4,4,0), N=30)
b = D.distrib1(a, 0.01)

test.testA([b], 1)
