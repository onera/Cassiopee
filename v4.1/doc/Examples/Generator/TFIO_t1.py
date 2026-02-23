# - TFIO (array) -
import Generator as G
import Geom as D
import KCore.test as test

a = D.circle((0,0,0), 1., N=41)
r = G.TFIO(a)
test.testA(r, 1)
