# - compressCartesian (pyTree) -
import Compressor.PyTree as Compressor
import Generator.PyTree as G
import KCore.test as test

# 3D
a = G.cart((0,0,0), (1,1,1), (10,10,10))
Compressor._compressCartesian(a)
test.testO(a, 1)

# 2D
a = G.cart((0,0,0), (1,1,1), (10,10,1))
Compressor._compressCartesian(a)
test.testO(a, 2)
