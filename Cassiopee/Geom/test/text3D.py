# - text3D (array) -
import Geom as D
import Converter as C

a = D.text3D("Cassiopee", smooth=1, thickness=2.)
C.convertArrays2File([a], 'out.plt')
