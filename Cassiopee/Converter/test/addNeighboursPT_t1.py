# - addNeighbours (pyTree) -
import Converter.PyTree as C
import Connector.PyTree as X
import Generator.PyTree as G
import Converter.elsAProfile as elsAProfile
import KCore.test as test

a = G.cylinder((0,0,0), 1., 1.5, 360., 0., 1., (30,30,10))
a = C.addBC2Zone(a, 'overlap1', 'BCOverlap', 'jmax')
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
b = G.cart((-10.,-10.,-10.),(0.4,0.4,0.4), (50,50,50))

t = C.newPyTree(['Cyl',a,'Cart',b])
t = X.applyBCOverlaps(t)
tp = elsAProfile.buildBCOverlap(t)
elsAProfile._rmGCOverlap(tp)
tp = elsAProfile.addNeighbours__(tp)
test.testT(tp, 1)
