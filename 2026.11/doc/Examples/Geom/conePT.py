# - cone (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C

a = D.cone((0,0,0), 1. , 0.5, 1.)
C.convertPyTree2File(a, 'out.cgns')
