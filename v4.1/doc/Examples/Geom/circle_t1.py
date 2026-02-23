# - circle (array) -
import Geom as D
import KCore.test as test

a = D.circle((0,0,0), 1. , 0., 360.)
test.testA([a], 1)
test.writeCoverage(100)
