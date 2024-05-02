# - sphere (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C

a = D.sphere((0,0,0), 1., 20)
C.convertPyTree2File(a, 'out.cgns')
