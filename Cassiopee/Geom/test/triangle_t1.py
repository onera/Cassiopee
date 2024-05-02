# - triangle -
import Geom as D
import KCore.test as test

a = D.triangle((0,0,0), (0.1,0.,0.1), (0.05, 0.08, 0.1))
test.testA([a],1)
test.writeCoverage(100)
