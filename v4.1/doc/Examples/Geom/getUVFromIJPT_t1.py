# - getUVFromIJ (pyTree) -
import Geom.PyTree as D
import Generator.PyTree as G
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (10,10,1))

D._getUVFromIJ(a)

test.testT(a, 1)
