# - airfoil (array) -
import Geom as D
import Converter as C

# list catalog
a = D.airfoil()

# list catalog of given series
a = D.airfoil("SIKORSKY")

# get one profile
a = D.airfoil("SIKORSKY/SIKORSKYSSC-A07AIRFOIL")
C.convertArrays2File(a, 'out.plt')
