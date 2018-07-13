# - F (array) -
import Converter as C
import Generator as G
import KCore.test as test

def F(x): return x

# Structure 1D
a = G.cart( (0,0,0), (1,1,1), (10,1,1) )
a = C.initVars(a, '{F}={x}+{y}')
b = F(a)
test.testA([b], 1)

# Structure 2D
a = G.cart( (0,0,0), (1,1,1), (10,10,1) )
a = C.initVars(a, '{F}={x}+{y}')
b = F(a)
test.testA([b], 2)

# Structure 3D
a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
a = C.initVars(a, '{F}={x}+{y}')
b = F(a)
test.testA([b], 3)

# BAR
a = G.cartTetra( (0,0,0), (1,1,1), (10,1,1) )
a = C.initVars(a, '{F}={x}+{y}')
b = F(a)
test.testA([b], 4)

# TRI
a = G.cartTetra( (0,0,0), (1,1,1), (10,10,1) )
a = C.initVars(a, '{F}={x}+{y}')
b = F(a)
test.testA([b], 5)

# QUAD
a = G.cartHexa( (0,0,0), (1,1,1), (10,10,1) )
a = C.initVars(a, '{F}={x}+{y}')
b = F(a)
test.testA([b], 6)

# TETRA
a = G.cartTetra( (0,0,0), (1,1,1), (10,10,10) )
a = C.initVars(a, '{F}={x}+{y}')
b = F(a)
test.testA([b], 7)

# HEXA
a = G.cartHexa( (0,0,0), (1,1,1), (10,10,10) )
a = C.initVars(a, '{F}={x}+{y}')
b = F(a)
test.testA([b], 8)

# PENTA
a = G.cartPenta( (0,0,0), (1,1,1), (10,10,10) )
a = C.initVars(a, '{F}={x}+{y}')
b = F(a)
test.testA([b], 9)

# PYRA

# NGON 1D
a = G.cartNGon( (0,0,0), (1,1,1), (10,1,1) )
a = C.initVars(a, '{F}={x}+{y}')
b = F(a)
test.testA([b], 11)

# NGON 2D
a = G.cartNGon( (0,0,0), (1,1,1), (10,10,1) )
a = C.initVars(a, '{F}={x}+{y}')
b = F(a)
test.testA([b], 12)

# NGON 3D
a = G.cartNGon( (0,0,0), (1,1,1), (10,10,10) )
a = C.initVars(a, '{F}={x}+{y}')
b = F(a)
test.testA([b], 13)

# liste d'arrays
a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
b = G.cartTetra( (0,0,0), (1,1,1), (10,10,10) )
A = F([a,b])
test.testA(A, 14)
