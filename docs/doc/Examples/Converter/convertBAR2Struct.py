# - convertBAR2Struct (array) -
import Converter as C
import Generator as G
import Geom as D
a = D.circle((0.,0.,0.),1.)
a = C.convertArray2Hexa(a); a = G.close(a)
b = C.convertBAR2Struct(a)
C.convertArrays2File([a,b], 'out.plt')
