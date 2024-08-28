# - cartr2 (array) -
import Generator as G
import KCore.test as test

a = G.cartr2((10,5,1), (1,1,1), (1.5,1.3,1.), (200.,100.,100.))
test.testA(a, 1)
