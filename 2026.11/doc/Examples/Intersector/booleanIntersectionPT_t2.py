# - boolean intersection (pyTree) -
# BAR
import Intersector.PyTree as XOR
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D
import KCore.test as test

c1 = D.circle( (0,0,0), 1, N=100 )
c2 = D.circle( (0.2,0,0), 1, N=50 )

c1 = C.convertArray2Tetra(c1); c1 = G.close(c1)
c2 = C.convertArray2Tetra(c2); c2 = G.close(c2)

x = XOR.booleanIntersection(c1, c2, tol=0.)
test.testT(x)
