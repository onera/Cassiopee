# - line (array) -
import Geom as D
import KCore.test as test

a = D.line((0,0,0), (1,0,0), N=100)
test.testA([a], 1)

test.writeCoverage(100)
