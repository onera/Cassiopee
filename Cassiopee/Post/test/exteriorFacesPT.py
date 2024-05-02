# - exteriorFaces (pyTree)-
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G

a = G.cartTetra((0,0,0), (1,1,1), (4,4,6))
b = P.exteriorFaces(a)
C.convertPyTree2File(b, 'out.cgns')
