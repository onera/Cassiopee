# - subzone (pyTree) -
import Converter.PyTree as C
import Transform.PyTree as T
import Generator.PyTree as G

a = G.cart((0,0,0), (1,1,1), (10,20,1))
a = T.subzone(a, (3,3,1), (7,8,1))
C.convertPyTree2File(a, 'out.cgns')
