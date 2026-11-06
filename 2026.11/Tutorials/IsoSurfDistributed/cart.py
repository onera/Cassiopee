# Case
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T

a = G.cart((0,0,0), (1,1,1), (100,100,50) )
a = C.initVars(a, 'F = sin({CoordinateX}*0.1)*cos({CoordinateY}*0.1)*sin({CoordinateZ}*0.1)')
a = T.splitSize(a, N=3000)
C.convertPyTree2File(a, 'cart.cgns')
