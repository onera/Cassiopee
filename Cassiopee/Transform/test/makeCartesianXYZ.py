# - makeCartesian(array) -
import Generator as G
import Transform as T
import Converter as C

a = G.cart((0.,0.,0.),(1.,1.,1.),(11,12,13))
a = T.reorder(a, (3,2,-1))
a = T.makeCartesianXYZ(a)
C.convertArrays2File(a, 'out.plt')
