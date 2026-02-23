# - cylinder (pyTree) -
import Geom.PyTree as D
import KCore.test as test

a = D.cylinder((0,0,0), 1., 10.)
b = D.cylinder((3,0,0), 1., 5., N=20, ntype='QUAD')
test.testT(a+[b], 1)
