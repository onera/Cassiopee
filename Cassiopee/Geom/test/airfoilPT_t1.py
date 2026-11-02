# - profile (pyTree) -
import Geom.PyTree as D
import KCore.test as test

# get one profile
a = D.airfoil("SIKORSKY/SIKORSKYSSC-A07AIRFOIL")
test.testT(a, 1)
