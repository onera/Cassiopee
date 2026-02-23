# - translate (pyTree) -
import Transform.PyTree as T
import Generator.PyTree as G
import Converter.PyTree as C

a = G.cart((0,0,0), (1,1,1), (10,10,3))
T._translate(a, (10.,0.,0.))
C.convertPyTree2File(a, 'out.cgns')
