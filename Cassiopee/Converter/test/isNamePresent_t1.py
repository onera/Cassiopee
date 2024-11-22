# - isNamePresent (array) -
import Converter as C
import Generator as G
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (50,50,50))
C._initVars(a, 'F', 1.)

b = G.cart((0,0,0), (1,1,1), (50,50,50))
C._initVars(b, 'centers:G', 2.)

r0 = C.isNamePresent(a, 'F')
r1 = C.isNamePresent([a, b], 'F')
r2 = C.isNamePresent([a, b], 'K')
test.testO([r0,r1,r2], 1)
