# - convertLO2HO (array) -
import Converter as C
import Geom as D

a = D.triangle((0,0,0), (1,0,0), (1,1,0))
a = C.convertLO2HO(a, mode=0)
# save to add when ready
