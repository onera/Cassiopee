# - isFinite (array) -
import Generator as G
import Converter as C
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (10,10,10))
a = C.initVars(a, 'F', 1.)
test.testO(C.isFinite(a), 1)
test.testO(C.isFinite(a, var='x'), 2)
test.testO(C.isFinite(a, var='F'), 3)
