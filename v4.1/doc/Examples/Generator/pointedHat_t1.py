# - pointedHat -
import Geom as D
import Generator as G
import KCore.test as test

# i-array
c = D.circle( (0,0,0), 1., 360., 0., 100)
surf = G.pointedHat(c,(0.,0.,1.))
test.testA([surf],1)

# ij-array
a = G.cylinder((0,0,0), 1., 1.5, 360., 0., 1., (100,30,1))
surf = G.pointedHat(a,(0.,0.,1.))
test.testA([surf],2)
