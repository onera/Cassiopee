# - profile (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C
import KCore.test as test

# get one profile
a = D.profile("SIKORSKY/SIKORSKYSSC-A07AIRFOIL")
C.convertPyTree2File(a, 'out.cgns')
test.testT(a, 1)
