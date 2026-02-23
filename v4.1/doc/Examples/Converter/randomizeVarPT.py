# - randomizeVar (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
a = G.cart((0,0,0),(1,1,1),(11,11,1))
C._randomizeVar(a, 'CoordinateZ', 0.1, 1.)
C.convertPyTree2File(a, "out.cgns")
