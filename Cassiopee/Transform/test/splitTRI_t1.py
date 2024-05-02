# - splitTRI (array) -
import Transform as T
import Generator as G
import Converter as C
import KCore.test as test

a = G.cart((0,0,0),(1,1,1),(5,5,1))
a = C.convertArray2Tetra(a)
c = [[10,16,22], [2,8,9]]
d = T.splitTRI(a, c)
test.testA(d,1)
