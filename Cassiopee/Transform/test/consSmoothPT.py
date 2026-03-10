# - consSmooth (pyTree) -
import Transform.PyTree as T
import Converter.PyTree as C
import Geom.PyTree as D
import Converter.Internal as I

sweeps = 5
l1 = D.line((0.,0.,0.), (0.,1.,0.), N=5)
l2 = D.line((0.,1.,0.), (1.,1.,0.), N=5)
l3 = D.line((1.,1.,0.), (1.,0.,0.), N=5)
l4 = D.line((1.,0.,0.), (0.,0.,0.), N=5)

a = T.join([l1,l2,l3,l4]); I.setName(a, 'carreOriginal')
b = T.consSmooth(a,sweeps); I.setName(b, 'carreLisse')

t = C.newPyTree(['Base',a,b])
C.convertPyTree2File(t, "out.cgns")