# - translate (array) -
import Transform as T
import Generator as G
import KCore.test as test

# in place array1

# STRUCT - ARRAY1 - OK
a = G.cart((0,0,0), (1,1,1), (5,5,5))
T._translate(a, (1.,0,0))
test.testA(a, 1)

# TETRA - ARRAY1 - OK
a = G.cartTetra((0,0,0), (1,1,1), (5,5,5))
T._translate(a, (1.,0,0))
test.testA(a, 2)


# NGON - ARRAY1 - OK
a = G.cartNGon((0,0,0), (1,1,1), (5,5,5))
T._translate(a, (1.,0,0))
test.testA(a, 3)
