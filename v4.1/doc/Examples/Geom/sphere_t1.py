# - sphere (array) -
import Geom as D
import KCore.test as test

a = D.sphere((0,0,0), 1., 20)
test.testA([a],1)
test.writeCoverage(100)
