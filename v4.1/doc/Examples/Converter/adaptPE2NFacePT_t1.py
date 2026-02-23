# - adaptNFace2PE (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G
import KCore.test as test

# NGon (CGNSv3)
a = G.cartNGon((0,0,0), (1,1,1), (10,10,10))
# a = C.initVars(a, '{G}={CoordinateX}')
# a = C.initVars(a, '{centers:F}={centers:CoordinateY}')
a2 = Internal.adaptNFace2PE(a, remove=True)
a2 = Internal.adaptPE2NFace(a, remove=False)
test.testT(a2,1)
a2 = Internal.adaptPE2NFace(a, remove=True)
test.testT(a2,2)

# NGon (CGNSv4)
b = G.cartNGon((0,0,0), (1,1,1), (10,10,10), api=3)
# b = C.initVars(b, '{G}={CoordinateX}')
# b = C.initVars(b, '{centers:F}={centers:CoordinateY}')
b2 = Internal.adaptNFace2PE(b, remove=True)
b2 = Internal.adaptPE2NFace(b, remove=False)
test.testT(b2,3)
b2 = Internal.adaptPE2NFace(b, remove=True)
test.testT(b2,3)
