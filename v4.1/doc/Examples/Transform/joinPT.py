# - join (pyTree) -
import Geom.PyTree as D
import Transform.PyTree as T
import Converter.PyTree as C

a1 = D.naca(12., 5001)
a2 = D.line((1.,0.,0.),(20.,0.,0.),5001)
a = T.join(a1, a2)
C.convertPyTree2File(a, "out.cgns")
