# - quadrangle (array) -
import Geom as D
import KCore.test as test

a = D.quadrangle((0,0,0), (0.1,0.,0.1), (0.05, 0.08, 0.1), (0.02,0.05,0.1))
test.testA([a],1)
test.writeCoverage(100)
