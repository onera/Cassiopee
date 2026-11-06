# - addkplane (pyTree) -
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.PyTree as C

a = G.cart((0.,0.,0.),(0.1,0.1,1.),(10,10,2))
a = T.addkplane(a)
C.convertPyTree2File(a, 'out.cgns')
