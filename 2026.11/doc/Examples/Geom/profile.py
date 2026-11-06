# - profile (array) -
import Geom as D
import Converter as C

# list catalog
a = D.profile()

# list catalog of given series
a = D.profile("SIKORSKY")

# get one profile
a = D.profile("SIKORSKY/SIKORSKYSSC-A07AIRFOIL")
C.convertArrays2File(a, 'out.plt')
