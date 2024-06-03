# - translate (pyTree) -
import Transform.PyTree as T
import Generator.PyTree as G
import Converter.PyTree as C
import KCore.test as test

# STRUCT - OK
a = G.cart((0,0,0), (1,1,1), (5,5,5))
T._translate(a, (1,0,0))
test.testT(a, 1)

# TETRA - OK
a = G.cartTetra((0,0,0), (1,1,1), (5,5,5))
T._translate(a, (1,0,0))
test.testT(a, 2)

# Multi ELEMENTS - OK
a1 = G.cartHexa((0,0,0), (1,1,1), (10,10,10))
a2 = G.cartHexa((9,0,0), (1,1,1), (10,10,10))
a = C.mergeConnectivity(a1, a2)
T._translate(a, (1,0,0))
test.testT(a, 3)

# NGONv3 - OK
a = G.cartNGon((0,0,0), (1,1,1), (5,5,5))
T._translate(a, (1,0,0))
test.testT(a, 4)

# NGONv4 - OK
a = G.cartNGon((0,0,0), (1,1,1), (5,5,5), api=3)
T._translate(a, (1,0,0))
test.testT(a, 5)

