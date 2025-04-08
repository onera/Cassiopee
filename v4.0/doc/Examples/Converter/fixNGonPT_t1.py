# - fixNGon (pyTree) -
import Generator.PyTree as G
import Converter.Internal as Internal
import KCore.test as test

# NGon avec PE sans NFace
a = G.cartNGon((0,0,0), (1,1,1), (10,10,10))
Internal._adaptNFace2PE(a)
Internal._fixNGon(a)
test.testT(a, 1)

# NGon avec PE decale du nombre de faces
a = G.cartNGon((0,0,0), (1,1,1), (10,10,10))
Internal._adaptNFace2PE(a, remove=False)
cFE = Internal.getNodeFromName2(a, 'ParentElements')
cFE[1] = cFE[1]+2430*(cFE[1]>0)
Internal._fixNGon(a)
test.testT(a, 2)

# NGonv4 (sort un NGonv3)
a = G.cartNGon((0,0,0), (1,1,1), (10,10,10), api=3)
Internal._fixNGon(a)
test.testT(a, 3)
