# - refine (pyTree) -
import Post.PyTree as P
import Converter.PyTree as C
import Generator.PyTree as G

# Linear with indicator field
a = G.cartTetra((0,0,0), (1,1,1), (10,10,1))
a = C.initVars(a, 'centers:indic', 1.)
a = P.refine(a, 'indic')
C.convertPyTree2File(a, 'out.cgns')
