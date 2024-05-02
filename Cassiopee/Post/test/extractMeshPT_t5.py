# - extractMesh (PyTree) -
# Tests pour mode=accurate
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G
import KCore.test as test

# extraction
b = G.cart((0.5,0.5,0.5),(1,1,1),(9,9,9))

# STRUCT
a = G.cart((0,0,0),(1,1,1),(10,10,10))
a = C.initVars(a, '{centers:Density}={centers:CoordinateX}')
b = P.extractMesh(a, b, mode='accurate')
test.testT(b, 1)

# TETRA
a = G.cartTetra((0,0,0),(1,1,1),(10,10,10))
a = C.initVars(a, '{centers:Density}={centers:CoordinateX}')
b = P.extractMesh(a, b, mode='accurate')
test.testT(b, 2)

# HEXA
a = G.cartHexa((0,0,0),(1,1,1),(10,10,10))
a = C.initVars(a, '{centers:Density}={centers:CoordinateX}')
b = P.extractMesh(a, b, mode='accurate')
test.testT(b, 2)

# NGON
a = G.cartNGon((0,0,0),(1,1,1),(10,10,10))
a = C.initVars(a, '{centers:Density}={centers:CoordinateX}')
b = P.extractMesh(a, b, mode='accurate')
test.testT(b, 2)
