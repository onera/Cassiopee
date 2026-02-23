# - convertArray2NGon(pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (2,2,1))
b = C.convertArray2NGon(a)
C.convertPyTree2File(b, 'out.cgns')
