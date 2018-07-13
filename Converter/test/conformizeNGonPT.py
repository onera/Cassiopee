# - conformizeNGon (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Transform.PyTree as T

a = G.cartNGon((0,0,0),(0.1,0.1,1),(11,11,1))
b = G.cartNGon((1.,0,0),(0.1,0.2,1),(11,6,1))
res = T.join(a,b)
res = C.conformizeNGon(res)
C.convertPyTree2File(res, 'out.cgns')
