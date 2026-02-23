# - splitCurvatureRadius (pyTree)-
import Converter.PyTree as C
import Geom.PyTree as D
import Transform.PyTree as T

a = D.naca(12.5000)
zones = T.splitCurvatureRadius(a, 10.)
C.convertPyTree2File(zones+[a], 'out.cgns')
