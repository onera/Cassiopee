# - splitCurvatureAngle (pyTree) -
import Converter.PyTree as C
import Geom.PyTree as D
import Transform.PyTree as T

a = D.naca(12,101)
a2 = D.line((1,0,0), (2,0,0), 50)
a = T.join(a, a2)
a2 = D.line((2,0,0), (1,0,0), 50)
a = T.join(a, a2)
zones = T.splitCurvatureAngle(a, 20.)
C.convertPyTree2File(zones+[a], 'out.cgns')
