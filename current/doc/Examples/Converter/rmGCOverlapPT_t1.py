# - rmGCOverlaps (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.elsAProfile as elsAProfile
import KCore.test as test

a = G.cylinder((0,0,0), 1., 1.5, 360., 0., 1., (30,30,10))
a = C.addBC2Zone(a, 'overlap1', 'BCOverlap', 'jmin')
t = C.newPyTree(['Base',a])
tp = elsAProfile.rmGCOverlap(t)
test.testT(tp, 1)
elsAProfile._rmGCOverlap(t)
test.testT(t, 1)
