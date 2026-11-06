# - distrib2 (array) -
import Geom as D
import KCore.test as test

a = D.line((0,0,0), (4,4,0), N=30)
b = D.distrib2(a, 0.001, 0.1, algo=0)
test.testA([b], 1)
b = D.distrib2(a, 0.001, 0.1, algo=1)
test.testA([b], 2)
