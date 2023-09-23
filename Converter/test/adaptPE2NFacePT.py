# - adaptPE2NFace (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G

# NGon
a = G.cartNGon((0,0,0), (1,1,1), (10,10,10))
a = Internal.adaptNFace2PE(a, remove=True)

a = Internal.adaptPE2NFace(a, remove=True)
C.convertPyTree2File(a, 'out.cgns')
