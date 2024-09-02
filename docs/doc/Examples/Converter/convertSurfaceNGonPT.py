# - convertSurfaceNGon (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G

# type A : NGON=bars, NFACE=polygon
a = G.cartNGon((0,0,0), (1,1,1), (10,10,1), api=3)
Internal.printTree(a)

# type B : NGON=polygon, NFACE=NULL
b = C.convertSurfaceNGon(a)
Internal.printTree(b)

# Back to type A
c = C.convertSurfaceNGon(b)
Internal.printTree(c)
