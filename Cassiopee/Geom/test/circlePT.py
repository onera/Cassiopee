# - circle (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C

a = D.circle((0,0,0), 1. , 0., 360.)
C.convertPyTree2File(a, 'out.cgns')
