# - compressFields (pyTree) -
import Compressor.PyTree as Compressor
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal

a = G.cart((0,0,0), (1,1,1), (10,10,10))
C._initVars(a, '{F}={CoordinateX}')
C._initVars(a, '{centers:G}={centers:CoordinateY}')
Compressor._compressFields(a, tol=1.e-6, ctype=5)
Compressor._uncompressAll(a)
Internal.printTree(a)
C.convertPyTree2File(a, 'out.cgns')
