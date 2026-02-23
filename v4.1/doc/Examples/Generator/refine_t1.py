# - refine (array) -
import Generator as G
import Converter as C
import KCore.test as test

a = G.cart((0,0,0), (0.1,0.1,0.1), (20,20,1))
a = G.refine(a, 1.5, 1)
test.testA([a], 1)
# liste de zones
a = G.cart((0,0,0), (0.1,0.1,0.1), (11,11,1))
b = G.cart((0,1,0), (0.1,0.1,0.1), (11,11,1))
A = [a,b]; A = C.initVars(A,'F',1.)
# facteur de raffinement non entier
B = G.refine(A,1.5,dir=1)
test.testA(B,2)

# facteur de raffinement entier
B = G.refine(A,3,dir=1)
test.testA(B,3)

# 3 directions
a = G.cart((0,0,0), (0.1,0.1,0.1), (11,11,1))
a = C.initVars(a,'F',1.)
b = G.refine(a,3,dir=0)
test.testA([b],4)

# 3D
a = G.cart((0,0,0), (0.1,0.1,0.1), (11,11,11))
a = C.initVars(a,'F',1.)
b = G.refine(a,3,dir=0)
test.testA([b],5)
