# - compressElements (pyTree) -
import Compressor.PyTree as Compressor
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal
import KCore.test as test

a = G.cartHexa((0,0,0), (1,1,1), (25,23,24))
Compressor._compressElements(a)
Compressor._uncompressAll(a)
test.testT(a, 1)

a = G.cartNGon((0,0,0), (1,1,1), (25,23,24))
Internal.printTree(a)

Compressor._compressElements(a)
Internal.printTree(a)
#Internal.printTree(a)

Compressor._uncompressAll(a)
Internal.printTree(a)
test.testT(a, 2)

