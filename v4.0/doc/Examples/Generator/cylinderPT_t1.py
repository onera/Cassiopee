# - cylinder (pyTree) -
import Generator.PyTree as G
import KCore.test as test

a = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0., 10., (50,50,30))
test.testT(a, 1)
test.writeCoverage(100)
