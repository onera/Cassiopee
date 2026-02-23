# - convertBAR2Struct (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Geom.PyTree as D
a = D.circle((0.,0.,0.),1.)
a = C.convertArray2Hexa(a); a = G.close(a)
b = C.convertBAR2Struct(a)
C.convertPyTree2File(b,'out.cgns')
