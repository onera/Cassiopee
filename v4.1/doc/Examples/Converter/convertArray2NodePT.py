# - convertArray2Node (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
a = C.convertArray2Node(a)
C.convertPyTree2File(a, 'out.cgns')
