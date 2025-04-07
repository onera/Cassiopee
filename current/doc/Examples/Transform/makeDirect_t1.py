# - makeDirect (array) -
import Generator as G
import Transform as T
import KCore.test as test

a = G.cart((0.,0.,0.),(1.,1.,1.),(10,10,10))
a = T.reorder(a, (1,2,-3)) # indirect now
a = T.makeDirect(a)
test.testA([a], 1)
