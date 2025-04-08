# - sphere6 (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C

A = D.sphere6((0,0,0), 1., 20)
b = D.sphere6((3,0,0), 1.2, N=20, ntype='QUAD')
C.convertPyTree2File(A+[b], 'out.cgns')
