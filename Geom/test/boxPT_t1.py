# - box (pyTree) -
import Geom.PyTree as D
import KCore.test as test

a = D.box((0,0,0), (1,1,1))
b = D.box((2,0,0), (3,1,1), N=30, ntype='QUAD')
test.testT(a+[b], 1)
