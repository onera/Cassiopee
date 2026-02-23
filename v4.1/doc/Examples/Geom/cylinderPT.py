# - cylinder (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C

a = D.cylinder((0,0,0), 1., 10.)
b = D.cylinder((3,0,0), 1., 5., N=20, ntype='QUAD')
C.convertPyTree2File(a+[b], 'out.cgns')
