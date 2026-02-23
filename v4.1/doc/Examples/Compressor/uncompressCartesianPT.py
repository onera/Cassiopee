# - uncompressCartesian (pyTree) -
import Compressor.PyTree as Compressor
import Generator.PyTree as G
import Converter.Internal as Internal

a = G.cart((0,0,0), (1,1,1), (10,10,10))
Compressor._compressCartesian(a)

Compressor._uncompressCartesian(a)
Internal.printTree(a)
