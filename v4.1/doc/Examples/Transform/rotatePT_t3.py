# - rotate (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T
import KCore.test as test

# Struct
z = G.cart((0,0,0), (1,1,1), (10,10,10))
T._rotate(z, (0,0,0), (0,0,1), 32.)
test.testT(z, 1)

# BE
z = G.cartHexa((0,0,0), (1,1,1), (10,10,10))
T._rotate(z, (0,0,0), (0,0,1), 32.)
test.testT(z, 2)

# ME
a = G.cartHexa((0,0,0), (1,1,1), (10,10,10))
b = G.cartTetra((9,0,0), (1,1,1), (10,10,10))
z = C.mergeConnectivity(a, b, boundary=0)
T._rotate(z, (0,0,0), (0,0,1), 32.)
test.testT(z, 3)

# NGON3
z = G.cartNGon((0,0,0), (1,1,1), (10,10,10), api=1)
T._rotate(z, (0,0,0), (0,0,1), 32.)
test.testT(z, 4)

# NGON4
z = G.cartNGon((0,0,0), (1,1,1), (10,10,10), api=3)
T._rotate(z, (0,0,0), (0,0,1), 32.)
test.testT(z, 5)
