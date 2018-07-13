# - cone -
import Geom as D
import KCore.test as test

a = D.cone((0,0,0), 1. , 0.5, 1.)
test.testA([a], 1)
