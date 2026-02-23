# - convertHO2LO (array) -
import Converter as C
import Geom as D

a = D.triangle((0,0,0), (1,0,0), (1,1,0))
a = C.convertLO2HO(a, mode=0)
a = C.convertHO2LO(a, mode=0)
C.convertArrays2File(a, 'out.plt')
