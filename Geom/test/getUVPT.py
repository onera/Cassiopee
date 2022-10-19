# - getUV (pyTree) -
import Geom.PyTree as D
import Generator.PyTree as G
import Converter.PyTree as C
import Post.PyTree as P

a = G.cartTetra((0,0,0), (1,1,1), (10,10,10))
a = P.exteriorFaces(a)
a = C.initVars(a, 'u', 0.)
a = C.initVars(a, 'v', 0.)
ret = D.getUV(a, 2., 0.)

C.convertPyTree2File(ret[0], 'out.cgns')
C.convertPyTree2File(ret[1], 'color.png')
C.convertPyTree2File(ret[2], 'bump.png')
