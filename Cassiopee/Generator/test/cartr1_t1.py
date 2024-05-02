# - cartr1 (array) -
import Generator as G
import KCore.test as test

a = G.cartr1((0,0,0), (1.,1.,1.), (1.1,1.2,1.), (10,10,10))
test.testA(a, 1)