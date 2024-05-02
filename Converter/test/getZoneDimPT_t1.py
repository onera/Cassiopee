# - getZoneDim (pyTree) -
import Generator.PyTree as G
import Geom.PyTree as D
import Converter.PyTree as C
import Converter.Internal as Internal
import KCore.test as test

# Structured
a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
dim = Internal.getZoneDim(a)
test.testO(dim, 1)

a = G.cart( (0,0,0), (1,1,0), (10,10,1) )
dim = Internal.getZoneDim(a)
test.testO(dim, 11)

a = G.cart( (0,0,0), (1,0,0), (10,1,1) )
dim = Internal.getZoneDim(a)
test.testO(dim, 12)

# NGON
a = G.cartNGon( (0,0,0), (1,1,1), (10,10,10), api=3 )
dim = Internal.getZoneDim(a)
test.testO(dim, 4)

a = G.cartNGon( (0,0,0), (1,1,0), (10,10,1), api=3 )
dim = Internal.getZoneDim(a)
test.testO(dim, 41)

a = G.cartNGon( (0,0,0), (1,0,0), (10,1,1), api=3 )
dim = Internal.getZoneDim(a)
test.testO(dim, 42)

# BE 3D
a = G.cartTetra( (0,0,0), (1,1,1), (10,10,10) )
dim = Internal.getZoneDim(a)
test.testO(dim, 2)

a = G.cartPyra( (0,0,0), (1,1,1), (10,10,10) )
dim = Internal.getZoneDim(a)
test.testO(dim, 21)

a = G.cartPenta( (0,0,0), (1,1,1), (10,10,10) )
dim = Internal.getZoneDim(a)
test.testO(dim, 22)

a = G.cartHexa( (0,0,0), (1,1,1), (10,10,10) )
dim = Internal.getZoneDim(a)
test.testO(dim, 23)

# BE 2D
a = G.cart( (0,0,0), (1,1,0), (10,10,1) )
a = C.convertArray2Hexa(a)
dim = Internal.getZoneDim(a)
test.testO(dim, 24)

a = G.cart( (0,0,0), (1,1,0), (10,10,1) )
a = C.convertArray2Tetra(a)
dim = Internal.getZoneDim(a)
test.testO(dim, 25)

# BE 1D & 0D
a = D.circle((0,0,0), 1. , 0., 360.)
a = C.convertArray2Hexa(a)
dim = Internal.getZoneDim(a)
test.testO(dim, 26)

a = G.cart( (0,0,0), (1,1,0), (10,10,1) )
a = C.convertArray2Node(a)
dim = Internal.getZoneDim(a)
test.testO(dim, 27)

# ME 3D
a = G.cartTetra( (0,0,0), (1,1,1), (10,10,10) )
b = G.cartHexa( (9,0,0), (1,1,1), (10,10,10) )
a = C.mergeConnectivity(a, b, boundary=0)
dim = Internal.getZoneDim(a)
test.testO(dim, 3)

# ME 2D
a = G.cart( (0,0,0), (1,1,0), (10,10,1) )
a = C.convertArray2Hexa(a)
b = G.cart( (9,0,0), (1,1,0), (10,10,1) )
b = C.convertArray2Tetra(b)
a = C.mergeConnectivity(a, b, boundary=0)
dim = Internal.getZoneDim(a)
test.testO(dim, 31)
