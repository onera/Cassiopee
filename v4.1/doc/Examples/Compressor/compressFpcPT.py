# - compressFpc (pyTree) -
import Compressor.PyTree as Compressor
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal

a = G.cart((0,0,0), (1,1,1), (3,3,3))
C._initVars(a, '{centers:F}=1.')

node = Internal.getNodeFromName(a, 'F')
Compressor._packNode(node, ctype=5)
Compressor._unpackNode(node)
C.convertPyTree2File(a, 'out.cgns')
