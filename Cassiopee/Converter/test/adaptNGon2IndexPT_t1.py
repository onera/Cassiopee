# - adaptNGon2Index (pyTree) -
import Generator.PyTree as G
import Converter.Internal as Internal
import KCore.test as test

a = G.cartNGon((0,0,0), (1,1,1), (10,10,10))
Internal._adaptNFace2PE(a, remove=False)
Internal._adaptNGon2Index(a)
Internal._adaptNFace2Index(a)
test.testT(a, 1)
