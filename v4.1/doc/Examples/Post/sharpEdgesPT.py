# - sharpEdges (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P
import Transform.PyTree as T
a1 = G.cart((0.,0.,0.),(1.5,1.,1.),(2,2,1))
a2 = T.rotate(a1,(0.,0.,0.),(0.,1.,0.),100.);a2[0] = 'cart2'
res = P.sharpEdges([a1,a2],alphaRef=45.)
t = C.newPyTree(['Edge',1,'Base',2]); t[2][1][2] += res; t[2][2][2]+=[a1,a2]
C.convertPyTree2File(t,"out.cgns")
