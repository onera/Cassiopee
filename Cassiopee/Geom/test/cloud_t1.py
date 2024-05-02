# - cloud (array) -
import Geom as D
import KCore.test as test
import numpy as np

n = 100
x = np.arange(n)
y = x
z = np.zeros_like(x)

a = D.cloud((x,y,z))
test.testA([a], 1)
test.writeCoverage(100)
