# - text3D (array) -
import Geom as D
import KCore.test as test

# deja teste dans text2D
a = D.text3D("ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789", smooth=1, font='text1')
test.testA([a], 1)
test.writeCoverage(100)
