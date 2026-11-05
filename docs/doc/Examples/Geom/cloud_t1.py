# - cloud (array) -
import Geom as D
import KCore.test as test
import numpy

n = 100
x = numpy.arange(n)
y = x
z = numpy.zeros_like(x)

a = D.cloud((x,y,z))
test.testA([a], 1)
test.writeCoverage(100)
