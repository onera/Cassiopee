# - mergeConnectivity (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

a = G.cartHexa( (0,0,0), (1,1,1), (10,10,10) )
b = G.cartHexa( (0,0,0), (1,1,1), (10,10,1) )

# merge boundary connectivity
c = C.mergeConnectivity(a, b, boundary=1)
test.testT(c, 1)

# merge element connectivity
b = G.cartTetra( (0,0,9), (1,1,1), (10,10,10) )
c = C.mergeConnectivity(a, b, boundary=0)
test.testT(c, 2)
