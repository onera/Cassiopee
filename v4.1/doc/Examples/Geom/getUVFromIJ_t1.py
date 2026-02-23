# - getUVFromIJ (array) -
import Geom as D
import Generator as G
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (10,10,1))

a = D.getUVFromIJ(a)

test.testA(a, 1)
