# - consSmooth (array) -
import Transform as T
import Converter as C
import Geom as D

sweeps = 5
l1 = D.line((0.,0.,0.), (0.,1.,0.), N=5)
l2 = D.line((0.,1.,0.), (1.,1.,0.), N=5)
l3 = D.line((1.,1.,0.), (1.,0.,0.), N=5)
l4 = D.line((1.,0.,0.), (0.,0.,0.), N=5)

a = T.join([l1,l2,l3,l4])
b = T.consSmooth(a, sweeps)

C.convertArrays2File([a,b], "out.plt")
