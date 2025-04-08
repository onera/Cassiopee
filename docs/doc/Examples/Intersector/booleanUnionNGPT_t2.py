# - adapts a cells with respect to b points (PyTree) -
#
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.PyTree as C
import Converter.Internal as I
import Intersector.PyTree as XOR
import Connector.PyTree as XC
import KCore.test as test

m = G.cart((0,0,0), (1,1,1), (20,20,20))
#C.convertPyTree2File(m, 'm.cgns')
t = C.newPyTree(['Base',m])
# split into 2 part
t = T.splitNParts(t, 2, multigrid=0, dirs=[1,2,3])

# turn one part to TETRA mesh
zones = I.getZones(t)
p1 = zones[0]
p1 = C.convertArray2Tetra(p1, split='withBarycenters')
p1 = C.convertArray2NGon(p1)
# C.convertPyTree2File(p1, 'p1.cgns')

p2 = C.convertArray2NGon(zones[1])
# C.convertPyTree2File(p2, 'p2.cgns')

x = XOR.booleanUnion(p1,p2, multi_zone=True) #conformize the join
#C.convertPyTree2File(x, 'u.cgns')

# compute the join
# x = XC.connectMatch(x, tol=1.e-6, dim=3)

test.testT(x,1)
