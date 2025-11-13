# - profile (array) -
import Geom as D
import Converter as C

# write catalog
a = D.profile()

# get one profile
a = D.profile("SIKORSKY/SIKORSKYSSC-A07AIRFOIL")
C.convertArrays2File(a, 'out.plt')
