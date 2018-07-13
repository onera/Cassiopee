# - empty (numpy like) -
import KCore
import KCore.test as test

# Return un numpy aligne (fortran, double)

a = KCore.empty((1200,3), 64)
a[:] = 1.
test.testO(a, 1)

a = KCore.empty(1200, 64)
a[:] = 1.
test.testO(a, 2)

a = KCore.empty((1200), 64)
a[:] = 1.
test.testO(a, 3)
