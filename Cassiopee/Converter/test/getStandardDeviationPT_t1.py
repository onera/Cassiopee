# - getStandardDeviation (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

a = G.cart( (0,0,0), (1.,1.,1.), (11,1,1) )
stdval = C.getStandardDeviation(a, 'CoordinateX')
test.testO(stdval, 1)
