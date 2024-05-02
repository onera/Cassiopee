# - deform (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T

a = G.cart((0.,0.,0.),(1.,1.,1.),(10,10,10))
C._initVars(a,'dx', 10.)
C._initVars(a,'dy', 0)
C._initVars(a,'dz', 0)
b = T.deform(a)
C.convertPyTree2File(b, 'out.cgns')
