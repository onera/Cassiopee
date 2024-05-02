# - getCurvatureHeight(pyTree) -
import Converter.PyTree as C
import Geom.PyTree as D

a = D.polyline([(0.,0.,0.),(1.,1.,0.),(2.,0.,0.)])
a = D.getCurvatureHeight(a)
C.convertPyTree2File(a, 'out.cgns')
