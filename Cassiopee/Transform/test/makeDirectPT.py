# - makeDirect (pyTree) -
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.PyTree as C

a = G.cart((0.,0.,0.),(1.,1.,1.),(10,10,10))
a = T.reorder(a, (1,2,-3)) # indirect now
a = T.makeDirect(a)
C.convertPyTree2File(a, 'out.cgns')
