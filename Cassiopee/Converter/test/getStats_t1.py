# - getStats (array) -
import Converter as C
import Generator as G
import KCore.test as test

a = G.cart((0,0,0), (1.,1.,1.), (11,1,1))
stats = C.getStats(a, 'x')
test.testO(stats, 1)