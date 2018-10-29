# - makeCartesian(pyTree) -
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.PyTree as C

a = G.cart((0.,0.,0.),(1.,1.,1.),(11,12,13))
a = T.reorder(a, (3,2,-1))
a = T.makeCartesianXYZ(a)
C.convertPyTree2File(a, 'out.cgns')
