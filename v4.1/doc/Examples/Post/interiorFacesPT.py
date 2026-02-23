# - interiorFaces (pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G

a = G.cartTetra((0,0,0), (1,1.,1), (20,20,1))
b = P.interiorFaces(a)
t = C.newPyTree(['Base',1,b])
C.convertPyTree2File(t, 'out.cgns')
