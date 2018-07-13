# - adaptNFace2PE (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G
import KCore.test as test

a = G.cartNGon((0,0,0), (1,1,1), (10,10,10))
a = C.initVars(a, '{G}={CoordinateX}')
a = C.initVars(a, '{centers:F}={centers:CoordinateY}')
a2 = Internal.adaptNFace2PE(a, remove=False)
test.testT(a2,1)
a2 = Internal.adaptNFace2PE(a, remove=True)
test.testT(a2,2)
