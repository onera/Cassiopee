# - profile (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C

# list catalog
a = D.airfoil()

# list catalog of given series
a = D.airfoil("SIKORSKY")

# get one profile
a = D.airfoil("SIKORSKY/SIKORSKYSSC-A07AIRFOIL")
C.convertPyTree2File(a, 'out.cgns')
