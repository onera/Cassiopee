# - getNCells (array) -
import Converter as C
import Generator as G
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (10,10,10))
n = C.getNCells(a)
test.testO(n, 1)

a = G.cartTetra((0,0,0), (1,1,1), (10,10,10))
n = C.getNCells(a)
test.testO(n, 2)
