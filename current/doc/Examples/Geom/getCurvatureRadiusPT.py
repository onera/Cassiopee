# - getCurvatureRadius (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C

a = D.circle((0,0,0), 1, 10, 0, 10)
a = D.getCurvatureRadius(a)
C.convertPyTree2File(a, 'out.cgns')
