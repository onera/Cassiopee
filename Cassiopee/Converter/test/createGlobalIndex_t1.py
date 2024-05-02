# - createGlobalIndex (array) -
import Converter as C
import Generator as G
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (10,10,10))
C._createGlobalIndex(a)
test.testA(a, 1)

a = G.cartTetra((0,0,0), (1,1,1), (10,10,10))
C._createGlobalIndex(a)
test.testA(a, 2)

a = G.cartNGon((0,0,0), (1,1,1), (10,10,10))
C._createGlobalIndex(a)
test.testA(a, 3)
