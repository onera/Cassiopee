# - isoSurf (pyTree) -
import Post.PyTree as P
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

a = G.cartTetra( (-20,-20,-20), (0.5,0.5,0.5), (50,50,50))
a = C.initVars(a, '{field}={CoordinateX}*{CoordinateX}+{CoordinateY}*{CoordinateY}+{CoordinateZ}')

iso = P.isoSurf(a, 'field', value=5.)
t = C.newPyTree(['Base',2]); t[2][1][2] += iso
test.testT(t, 1)
