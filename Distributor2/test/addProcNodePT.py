# - addProcNode (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Distributor2.PyTree as D2

a = G.cart((0,0,0), (1,1,1), (10,10,10))
a = D2.addProcNode(a, 12)
C.convertPyTree2File(a, 'out.cgns')
