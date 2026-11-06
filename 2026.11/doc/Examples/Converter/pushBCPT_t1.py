# - pushBC (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import KCore.test as test

# STRUCT -> NGON
a = G.cart((0,0,0), (1,1,1), (10,10,10))
a = C.addBC2Zone(a, 'wall', 'BCWall', 'kmin')
b = G.cartNGon((0,0,0), (1,1,1), (10,10,10))
b = C.pushBC(a, b)
test.testT(b, 1)

# STRUCT -> BE (FACES)
b = G.cartHexa( (0,0,0), (1,1,1), (10,10,10) )
b = C.pushBC(a, b)
test.testT(b, 2)

# STRUCT -> BE (BCC)
b = G.cartHexa( (0,0,0), (1,1,1), (10,10,10) )
b = C.pushBC(a, b, type='BCC')
test.testT(b, 3)

# BE (BCC) -> NGON
a = G.cartHexa( (0,0,0), (1,1,1), (10,10,10) )
bc = G.cartHexa( (0,0,0), (1,1,1), (10,10,1) )
a = C.addBC2Zone(a, 'wall', 'BCWall', subzone=bc)
b = G.cartNGon( (0,0,0), (1,1,1), (10,10,10) )
b = C.pushBC(a, b)
test.testT(b, 4)

# BE (BCC) -> BE (F)
b = G.cartHexa( (0,0,0), (1,1,1), (10,10,10) )
b = C.pushBC(a, b)
test.testT(b, 5)

# NGON -> NGON
a = G.cartNGon( (0,0,0), (1,1,1), (10,10,10) )
bc = G.cart( (0,0,0), (1,1,1), (10,10,1) )
a = C.addBC2Zone(a, 'wall', 'BCWall', subzone=bc)

b = G.cartNGon( (0,0,0), (1,1,1), (10,10,10) )
b = C.pushBC(a, b)
test.testT(b, 6)
