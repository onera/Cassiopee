# - initConst (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Initiator.PyTree as I
import KCore.test as test

# STRUCT
z = G.cart((0.,0.,0.), (0.05,0.05,1.), (100,100,2))
I._initConst(z, MInf=0.8, loc='centers')
test.testT(z, 1)

# NGONv3
z = G.cartNGon((0,0,0), (0.05,0.05,1), (100,100,2), api=1)
I._initConst(z, MInf=0.8, loc='centers')
test.testT(z, 2)

# NGONv4
z = G.cartNGon((0,0,0), (0.05,0.05,1), (100,100,2), api=3)
I._initConst(z, MInf=0.8, loc='centers')
test.testT(z, 3)

# HEXA
z = G.cartHexa((0,0,0), (0.05,0.05,1), (100,100,2))
I._initConst(z, MInf=0.8, loc='centers')
test.testT(z, 4)

# ME
a = G.cartHexa((0,0,0), (0.05,0.05,1), (50,100,2))
b = G.cartTetra((2.45,0,0), (0.05,0.05,1), (50,100,2))
z = C.mergeConnectivity(a, b, boundary=0)
I._initConst(z, MInf=0.8, loc='centers')
test.testT(z, 5)
