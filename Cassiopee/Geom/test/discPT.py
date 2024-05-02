# - disc (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C

a = D.disc((0,0,0), 1.)
b = D.disc((3,0,0), 1., N=20, ntype='QUAD')
C.convertPyTree2File(a+[b], 'out.cgns')
