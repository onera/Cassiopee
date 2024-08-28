# - torus (array) -
import Geom as D
import KCore.test as test

# half torus
a = D.torus((0,0,0), 5., 2., 0., 180.,0.,360.,100, 50)
test.testA([a],1)

# trimmed torus
a = D.torus((0,0,0), 5., 2., 0., 120., 0., 90., 100, 36)
test.testA([a],2)
test.writeCoverage(100)
