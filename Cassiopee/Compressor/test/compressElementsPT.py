# - compressElements (pyTree) -
import Compressor.PyTree as Compressor
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal

a = G.cartHexa((0,0,0), (1,1,1), (25,23,24))
Compressor._compressElements(a)
Internal.printTree(a)
Compressor._uncompressAll(a)
Internal.printTree(a)
C.convertPyTree2File(a, 'out.cgns')
