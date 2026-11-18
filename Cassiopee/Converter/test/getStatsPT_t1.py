# - getStats (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

a = G.cart( (0,0,0), (1.,1.,1.), (11,1,1) )
stats = C.getStats(a, 'CoordinateX')
test.testO(stats, 1)
