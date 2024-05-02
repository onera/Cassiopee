# - adapt2FastP (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal

a = G.cartNGon((0,0,0), (1,1,1), (30,30,20))
Internal._adapt2FastP(a)
C.convertPyTree2File(a, 'out.cgns')
