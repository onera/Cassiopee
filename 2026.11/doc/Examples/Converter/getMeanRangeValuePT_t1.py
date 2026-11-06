# - getMeanRangeValue (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

a = G.cart((0,0,0), (1.,1.,1.), (11,1,1))
# get the mean of the 30% smallest values
meanval = C.getMeanRangeValue(a, 'CoordinateX', 0., 0.3)
test.testO(meanval, 1)
