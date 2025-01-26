# - convertArray2Tetra -
import Converter as C
import Generator as G
import KCore.test as test

a = G.cartPenta((0.,0.,0.), (1,1,1), (10,10,3))
a = C.convertArray2Tetra(a)
test.testA([a], 1)
