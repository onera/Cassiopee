# - getAngleRegularityMap (array) -
import Generator as G
import Converter as C
import Transform as T

a = G.cart((0,0,0), (1,1,1), (50,50,1))
a = T.deformPoint(a, (25,25,0), (1.,1.,0.), 2., 2.)
ac = C.node2Center(a)
reg = G.getAngleRegularityMap(a)
reg = C.addVars([ac,  reg])
C.convertArrays2File([reg], "out.plt")
