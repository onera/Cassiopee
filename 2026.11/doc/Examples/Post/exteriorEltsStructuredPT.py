# - exteriorEltsStructured (pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G

a = G.cart((0,0,0), (1,1,1), (10,10,1))
C._initVars(a, '{F}={CoordinateX}')
p = P.exteriorEltsStructured(a, 2)
C.convertPyTree2File(p, 'out.cgns')
