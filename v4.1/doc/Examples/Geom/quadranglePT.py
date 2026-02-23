# - quadrangle (PyTree) -
import Geom.PyTree as D
import Converter.PyTree as C

a = D.quadrangle((0,0,0.1), (0.1,0.,0.1), (0.05, 0.08, 0.1), (0.02,0.05,0.1))
b = D.quadrangle((0,0,0.1), (0.1,0.,0.1), (0.05, 0.08, 0.1), (0.02,0.05,0.1), N=20, ntype='QUAD')
C.convertPyTree2File(a+[b], 'out.cgns')
