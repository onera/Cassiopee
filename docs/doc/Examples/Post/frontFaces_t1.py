# - frontFaces (array) -
import Converter as C
import Generator as G
import Post as P
import KCore.test as test

def F(x, y, z):
    if (x + 2*y + z   > 20.): return 1
    else: return 0

# STRUCT
a = G.cart( (0,0,0), (1,1,1), (11,11,11) )
a = C.initVars(a, 'tag', F, ['x', 'y', 'z'])
t = C.extractVars(a, ['tag'])
f = P.frontFaces(a, t)
test.testA([f], 1)

# HEXA
a = G.cartHexa( (0,0,0), (1,1,1), (11,11,11) )
a = C.initVars(a, 'tag', F, ['x', 'y', 'z'])
t = C.extractVars(a, ['tag'])
f = P.frontFaces(a, t)
test.testA([f], 2)

# TETRA
a = G.cartTetra( (0,0,0), (1,1,1), (11,11,11) )
a = C.initVars(a, 'tag', F, ['x', 'y', 'z'])
t = C.extractVars(a, ['tag'])
f = P.frontFaces(a, t)
test.testA([f], 3)

# TRI
a = G.cartTetra( (0,0,0), (1,1,1), (11,11,1) )
a = C.initVars(a, 'tag', F, ['x', 'y', 'z'])
t = C.extractVars(a, ['tag'])
f = P.frontFaces(a, t)
test.testA([f], 4)

# QUAD
a = G.cartHexa( (0,0,0), (1,1,1), (11,11,1) )
a = C.initVars(a, 'tag', F, ['x', 'y', 'z'])
t = C.extractVars(a, ['tag'])
f = P.frontFaces(a, t)
test.testA([f], 5)
