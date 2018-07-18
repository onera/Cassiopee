# - point (array) -
import Geom as D
import KCore.test as test

a = D.point((0,0,0))
test.testA([a], 1)
test.writeCoverage(100)
