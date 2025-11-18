# - profile (pyTree) -
import Geom.PyTree as D
import KCore.test as test

# get one profile
a = D.profile("SIKORSKY/SIKORSKYSSC-A07AIRFOIL")
print(a)
test.testT(a, 1)
