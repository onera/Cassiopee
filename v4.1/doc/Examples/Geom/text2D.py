# - text2D (array) -
import Geom as D
import Converter as C

a = D.text2D("ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789", smooth=0, offset=1.)
C.convertArrays2File([a], 'out.plt')
