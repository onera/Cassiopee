# - getVolume (pyTree) -
import Geom.PyTree as D
import KCore.test as test

a = D.sphere((0,0,0), 1., N=30)
vol = D.getVolume(a)
test.testO(vol, 1)
