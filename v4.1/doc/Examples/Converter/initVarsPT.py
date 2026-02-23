# - initVars (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cart((0,0,0), (1,1,1), (10,10,10))
a = C.initVars(a, 'F', 0.)
a = C.initVars(a, 'centers:G', 1.)
C.convertPyTree2File(a, 'out.cgns')
