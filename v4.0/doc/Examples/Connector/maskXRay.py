# - maskXRay (array) -
import Converter as C
import Connector as X
import Geom as D
import Transform as T

surf = D.sphere((0,0,0), 0.5, 20)
surf = T.rotate(surf,(0.,0.,0.),(0.,1.,0.), 90.)
res = X.maskXRay__([surf])
C.convertArrays2File([res], "out.plt")
