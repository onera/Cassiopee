# - getUV (pyTree) -
import Geom.PyTree as D
import Generator.PyTree as G
import Converter.PyTree as C
import Post.PyTree as P

a = G.cartTetra((0,0,0), (1,1,1), (10,10,10))
a = P.exteriorFaces(a)
(a, color, normal) = D.getUV(a, 2., 1920)

C.convertPyTree2File(a, 'out.cgns')
C.convertPyTree2File(color, 'color.png')
C.convertPyTree2File(normal, 'bump.png')
