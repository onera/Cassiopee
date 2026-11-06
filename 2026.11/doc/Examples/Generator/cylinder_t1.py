# - cylinder (array) -
import Generator as G
import KCore.test as test

# 3D
a = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0., 10., (50,50,30))
test.testA([a], 1)

# 2D
a = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0., 10., (50,50,1))
test.testA([a], 2)

# 1D
a = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0., 10., (50,1,1))
test.testA([a], 3)
