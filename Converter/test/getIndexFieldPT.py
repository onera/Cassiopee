# - getIndexField (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cart((0,0,0), (1,1,1), (10,1,1))
C._getIndexField(a)
C.convertPyTree2File(a, 'out.cgns')
