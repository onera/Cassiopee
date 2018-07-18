# test pierce points XRAY (array)
# cas surface 2D avec body en BAR
import Converter as C
import Connector as X
import Geom as D

surf = D.circle((0,0,0), 0.5, 0., 360.)
surf = C.convertArray2Tetra(surf)
res = [surf]
res0 =  X.maskXRay__(res, 0.,2)
C.convertArrays2File([res0], "out.plt")
