# - profile (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C
import KCore.test as test

LOCAL = test.getLocal()

# get one profile
a = D.profile("SIKORSKY/SIKORSKYSSC-A07AIRFOIL")
#C.convertPyTree2File(a, LOCAL+'/out.cgns')
test.testT(a, 1)
