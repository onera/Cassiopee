# - isNamePresent (PyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (50,50,50))
C._initVars(a, 'F', 1.)
C._initVars(a, 'centers:G', 0.)

b = G.cart((0,0,0), (1,1,1), (50,50,50))
C._initVars(b, 'F', 2.)
C._initVars(b, 'centers:H', 3.)

t = C.newPyTree(['Base',a,b])

r0 = C.isNamePresent(a, 'F')
r1 = C.isNamePresent(a, 'centers:F')
r2 = C.isNamePresent([a, b], 'F')
r3 = C.isNamePresent([a, b], 'centers:G')
test.testO([r0,r1,r2,r3], 1)
