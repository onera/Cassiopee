# - convertSurfaceNGon (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G
import KCore.test as test

## NGon v3 ##############################

# type A : NGON=bars, NFACE=polygon
a = G.cartNGon((0,0,0), (1,1,1), (10,10,1), api=2)

# type B : NGON=polygon, NFACE=NULL
b = C.convertSurfaceNGon(a)
#test.testT(b, 1)

# type A : NGON=bars, NFACE=polygon
c = C.convertSurfaceNGon(b)
test.testT(c, 2)

## NGon v4 ##############################

# type A : NGON=bars, NFACE=polygon
a = G.cartNGon((0,0,0), (1,1,1), (10,10,1), api=3)

# type B : NGON=polygon, NFACE=NULL
b = C.convertSurfaceNGon(a)
#test.testT(b, 3)

# type A : NGON=bars, NFACE=polygon
c = C.convertSurfaceNGon(b)
test.testT(c, 4)
