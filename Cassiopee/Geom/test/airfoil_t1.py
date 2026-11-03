# - profile (array) -
import Geom as D
import KCore.test as test

# get one profile
a = D.airfoil("SIKORSKY/SIKORSKYSSC-A07AIRFOIL")
test.testA(a, 1)
