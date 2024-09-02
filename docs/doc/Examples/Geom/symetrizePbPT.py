# - symetrizePb (pyTree) -
import Converter.PyTree as C
import Geom.IBM as D_IBM
import Geom.PyTree as D
# Geometry
a = D.cone((0,0,0), 1. , 0., 1.)
D_IBM._setDfar(a, 10.)
D_IBM._setSnear(a,0.01)
tb = C.newPyTree(["BODY",a])
D_IBM._symetrizePb(tb,"BODY", 0.1, dir_sym=3)
C.convertPyTree2File(tb, 'out.cgns')
