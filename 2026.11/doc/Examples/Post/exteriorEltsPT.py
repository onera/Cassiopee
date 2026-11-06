# - exteriorElts (pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G

a = G.cartTetra((0,0,0), (1,1,1), (10,10,10))
b = P.exteriorElts(a)
C.convertPyTree2File(b, 'out.cgns')
