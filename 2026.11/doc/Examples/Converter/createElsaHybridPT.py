# - createElsaHybrid (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.elsAProfile as elsAProfile

# Must be NGon
a = G.cartNGon((0,0,0), (1,1,1), (10,10,10))
a = C.fillEmptyBCWith(a, 'farfield', 'BCFarfield')
elsAProfile._createElsaHybrid(a)
C.convertPyTree2File(a, 'out.cgns')
