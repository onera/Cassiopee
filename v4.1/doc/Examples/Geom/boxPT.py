# - box (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C

a = D.box((0,0,0), (1,1,1))
b = D.box((2,0,0), (3,1,1), N=30, ntype='QUAD')
C.convertPyTree2File(a+[b], 'out.cgns')
