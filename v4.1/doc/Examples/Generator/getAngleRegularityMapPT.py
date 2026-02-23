# - getAngleRegularityMap (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Transform.PyTree as T

a = G.cart((0,0,0), (1,1,1), (50,50,1))
a = T.deformPoint(a, (25,25,0), (1.,1.,0.), 2., 2.)
a = G.getAngleRegularityMap(a)
C.convertPyTree2File(a, 'out.cgns')
