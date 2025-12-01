# - profile (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C

# list catalog
a = D.profile()

# list catalog of given series
a = D.profile("SIKORSKY")

# get one profile
a = D.profile("SIKORSKY/SIKORSKYSSC-A07AIRFOIL")
C.convertPyTree2File(a, 'out.cgns')
