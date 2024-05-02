# - triangle (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C

a = D.triangle((0,0,0), (0.1,0.,0.1), (0.05, 0.08, 0.1))
C.convertPyTree2File(a, 'out.cgns')
