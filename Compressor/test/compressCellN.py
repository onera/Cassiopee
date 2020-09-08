# - compressCellN (pyTree) -
import Compressor.PyTree as Compressor
import Compressor.compressor as Co
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal
import numpy

a = G.cart((0,0,0), (1,1,1), (10,11,12))
C._initVars(a, '{centers:cellN}=1.')
Compressor._compressCellN(a)
Compressor._uncompressAll(a)
Internal.printTree(a)
C.convertPyTree2File(a, 'out.cgns')
