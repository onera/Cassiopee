# - selectCells2 (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P

a = G.cart((0,0,0), (1,1,1), (11,11,11))
a = C.initVars(a, 'tag', 1.)
b = P.selectCells2(a, 'tag')
C.convertPyTree2File(b, 'out.cgns')
