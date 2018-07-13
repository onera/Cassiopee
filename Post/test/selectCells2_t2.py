# - selectCells2 (array) -
import Converter as C
import Generator as G
import Post as P
import KCore.test as test

def F(x, y, z):
    if (x+2*y+z > 20.): return True
    else: return False

# Tetra
a = G.cartTetra( (0,0,0), (1,1,1), (11,11,11) )
tag = C.initVars(a, 'tag', F, ['x','y','z'])
tag = C.extractVars(tag,['tag'])
a = P.selectCells2(a, tag)
test.testA([a], 1)

# Penta
a = G.cartPenta( (0,0,0), (1,1,1), (11,11,11) )
tag = C.initVars(a, 'tag', F, ['x','y','z'])
tag = C.extractVars(tag,['tag'])
a = P.selectCells2(a, tag)
test.testA([a], 2)

# Pyra
a = G.cartPyra( (0,0,0), (1,1,1), (11,11,11) )
tag = C.initVars(a, 'tag', F, ['x','y','z'])
tag = C.extractVars(tag,['tag'])
a = P.selectCells2(a, tag)
test.testA([a], 3)

# NGon
a = G.cartNGon( (0,0,0), (1,1,1), (11,11,11) )
tag = C.initVars(a, 'tag', F, ['x','y','z'])
tag = C.extractVars(tag,['tag'])
a = P.selectCells2(a, tag)
test.testA([a], 4)

# Hexa
a = G.cartHexa( (0,0,0), (1,1,1), (11,11,11) )
tag = C.initVars(a, 'tag', F, ['x','y','z'])
tag = C.extractVars(tag,['tag'])
a = P.selectCells2(a, tag)
test.testA([a], 5)

# Quad
a = G.cartHexa( (0,0,0), (1,1,1), (11,11,1) )
tag = C.extractVars(tag,['tag'])
tag = C.initVars(a, 'tag', F, ['x','y','z'])
tag = C.extractVars(tag,['tag'])
a = P.selectCells2(a, tag)
test.testA([a], 6)

# Tri
a = G.cartTetra( (0,0,0), (1,1,1), (11,11,1) )
tag = C.initVars(a, 'tag', F, ['x','y','z'])
tag = C.extractVars(tag,['tag'])
a = P.selectCells2(a, tag)
test.testA([a], 7)

# Bar
a = G.cartTetra( (0,9,0), (1,1,1), (11,1,1) )
tag = C.initVars(a, 'tag', F, ['x','y','z'])
tag = C.extractVars(tag,['tag'])
a = P.selectCells2(a, tag)
test.testA([a], 8)
