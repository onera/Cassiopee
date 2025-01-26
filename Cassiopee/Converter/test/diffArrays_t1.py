# - diffArrays -
import Converter as C
import Generator as G
import KCore.test as T

ni = 11; nj = 11; nk = 11

# Test avec des variables differentes
a = G.cart( (0,0,0), (1,1,1), (ni,nj,nk) )
a = C.addVars(a, "F")
a = C.initVars(a, "F", 1.)
a = C.addVars(a, "Q")
a = C.initVars(a, "Q", 1.01)

b = G.cart( (0,0,0), (1,1,1), (ni,nj,nk) )
b = C.addVars(b, "F")
b = C.initVars(b, "F", 2.)
b = C.addVars(b, "Q")
b = C.initVars(b, "Q", 1.03)
b = C.addVars(b, "R")
b = C.initVars(b, "R", 1.12)

ret = C.diffArrays([a], [b])
T.testA(ret, 1)

# Test avec des coordonnees differentes
a = G.cart( (0.2,0.3,0), (1,1,1), (ni,nj,nk) )
a = C.addVars(a, "F")
a = C.initVars(a, "F", 1.)
a = C.addVars(a, "Q")
a = C.initVars(a, "Q", 1.01)

b = G.cart( (0,0,0), (1,1,1), (ni,nj,nk) )
b = C.addVars(b, "F")
b = C.initVars(b, "F", 2.)
b = C.addVars(b, "Q")
b = C.initVars(b, "Q", 1.03)

ret = C.diffArrays([a], [b])
T.testA(ret, 2)

# Test sans coordonnees
a = C.array('F,G',ni,nj,nk)
b = C.array('F,G',ni,nj,nk)
b = C.initVars(b, "F", 1.2)

ret = C.diffArrays([a], [b])
T.testA(ret, 3)
