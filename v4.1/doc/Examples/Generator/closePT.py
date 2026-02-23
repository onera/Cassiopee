# - close (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a1 = G.cart((0,0,0), (1,1,1), (10,10,1))
a2 = G.cart((9+1.e-2,0,0), (1,1,1), (10,10,1))
a3 = G.cart((0,-5.01,0),(1,1,1),(19,6,1))
a4 = G.cart((0,9.0001,0),(1,1,1),(10,6,1))
a5 = G.cart((9.01,9.0002,0),(1,1,1),(10,6,1))
t = C.newPyTree(['Base',2,a1,a2,a3,a4,a5])
t = G.close(t, 1.e-1)
C.convertPyTree2File(t, 'out.cgns')
