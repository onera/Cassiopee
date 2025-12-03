# - profile (array) -
import Geom as D
import KCore.test as test

# get one profile
a = D.profile("SIKORSKY/SIKORSKYSSC-A07AIRFOIL")
test.testA(a, 1)
