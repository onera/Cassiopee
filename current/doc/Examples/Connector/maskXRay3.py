# - maskXRay : points de percage (array) -
# cas 2D TRI

import Converter as C
# test pierce points XRAY
# cas surface 2D avec body en TRI

import Connector as X
import Generator as G
import Transform as T

surf = G.cylinder((0.,0.,0.), 0., 1., 360., 0., 1., (50,50,2))
surf = T.subzone(surf,(1,50,1),(50,50,2))
surf = C.convertArray2Tetra(surf)
res = [surf]
res0 =  X.maskXRay__(res, 0.,2)
C.convertArrays2File([res0], "out.plt")
