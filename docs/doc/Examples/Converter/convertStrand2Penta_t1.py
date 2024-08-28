# - convertSTrand2Penta (array) -
import Converter as C
import Generator as G
import KCore.test as test

a = G.cartPenta((0,0,0), (1,1,1), (3,3,3))

b = C.convertPenta2Strand(a)

c = C.convertStrand2Penta(b)

test.testA(c, 1)
