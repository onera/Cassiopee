# - getFamilyZones (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

a = G.cylinder((0,0,0), 1., 1.5, 0., 360., 1., (80,30,2))
b = G.cylinder((3,0,0), 1., 1.5, 0., 360., 1., (80,30,2))
c = G.cart((-0.1,0.9,0), (0.01,0.01,1.), (20,20,2))

C._tagWithFamily(a, 'CYLINDER')
C._tagWithFamily(b, 'CYLINDER')
C._tagWithFamily(c, 'CART')

t = C.newPyTree(['Base',a,b,c])
C._addFamily2Base(t[2][1], 'CYLINDER')
C._addFamily2Base(t[2][1], 'CART')
zones = C.getFamilyZones(t, 'CYLINDER')
test.testT(zones, 1)
