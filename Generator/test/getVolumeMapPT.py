# - getVolumeMap (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C

a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,3))
a = G.getVolumeMap(a)
C.convertPyTree2File(a, 'out.cgns')

# method=1
a = C.convertArray2NGon(a)
a = G.getVolumeMap(a, method=1, tol=1e-12)
C.convertPyTree2File(a, 'out2.cgns')
