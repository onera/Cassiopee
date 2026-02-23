# - uniformize (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C

a = D.polyline([(0,0,0), (1,1,0), (2,0,0), (3,1,0), (4,0,0)])
a = D.uniformize(a, N=100)

C.convertPyTree2File(a, 'out.cgns')
