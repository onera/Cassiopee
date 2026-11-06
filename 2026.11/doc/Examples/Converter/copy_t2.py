# - copy (array2) -
import Converter as C
import Generator as G
import KCore.test as test

# STRUCT - ARRAY2 - OK
a = G.cart((0,0,0), (1,1,1), (5,5,5), api=2)
b = C.copy(a)
test.testA(b, 1)

# TETRA - ARRAY2 - OK
a = G.cartTetra((0,0,0), (1,1,1), (5,5,5), api=2)
b = C.copy(a)
test.testA(b, 2)

# NGON - ARRAY2 - OK
a = G.cartNGon((0,0,0), (1,1,1), (5,5,5), api=2)
b = C.copy(a)
test.testA(b, 3)