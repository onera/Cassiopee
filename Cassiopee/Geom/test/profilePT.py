# - profile (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C

# write catalog
a = D.profile()

# get one profile
a = D.profile("SIKORSKY/SIKORSKYSSC-A07AIRFOIL")
C.convertPyTree2File(a, 'out.cgns')
