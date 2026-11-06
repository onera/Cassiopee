# - enforceX, enforceY, ... monotonic (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C

a = G.cart((0,0,0), (1,1,1), (20,20,10))
b = G.enforceX(a, 5., 0.2, 10, 5 )
b = G.enforceX(b, 15., 0.2, 10, 5 )
b = G.enforceY(b, 10., 0.1, 5, 5 )
C.convertPyTree2File(b, 'out.cgns')
