# - sliceNGonFaces (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

# - NGons -
a = G.cartNGon((2,0,0), (0.1,0.1,1), (10,10,2))
faceIndices = [1, 11, 21, 31, 41]
faceVertices, faceOffset = C.sliceNGonFaces(a, indices=faceIndices)
test.testO(faceVertices, 1)
test.testO(faceOffset, 2)
