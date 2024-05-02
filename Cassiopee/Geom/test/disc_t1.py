# - disc (array) -
import Geom as D
import KCore.test as test

a = D.disc((0,0,0), 1.)
test.testA(a, 1)
