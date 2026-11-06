# - boolean difference (array) -
# BAR
import Intersector as XOR
import Generator as G
import Converter as C
import Geom as D
import KCore.test as test

c1 = D.circle( (0,0,0), 1, N=100 )
c2 = D.circle( (0.2,0,0), 1, N=50 )

c1 = C.convertArray2Tetra(c1); c1 = G.close(c1)
c2 = C.convertArray2Tetra(c2); c2 = G.close(c2)

x = XOR.booleanMinus(c1, c2, tol=0.)
test.testA([x],1)
