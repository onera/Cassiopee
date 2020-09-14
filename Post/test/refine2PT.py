# - refine (pyTree) -
import Post.PyTree as P
import Converter.PyTree as C
import Generator.PyTree as G

# Refine with butterfly interpolation
a = G.cartTetra((0,0,0), (1,1,1), (10,10,1))
a = P.refine(a, w=1./64.)
C.convertPyTree2File(a, 'out.cgns')
