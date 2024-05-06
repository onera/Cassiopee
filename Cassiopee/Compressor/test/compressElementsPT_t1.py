# - compressElements (pyTree) -
import Compressor.PyTree as Compressor
import Generator.PyTree as G
import KCore.test as test

a = G.cartHexa((0,0,0), (1,1,1), (25,23,24))
Compressor._compressElements(a)
Compressor._uncompressAll(a)
test.testT(a, 1)

a = G.cartNGon((0,0,0), (1,1,1), (25,23,24))
Compressor._compressElements(a)
Compressor._uncompressAll(a)
test.testT(a, 2)
