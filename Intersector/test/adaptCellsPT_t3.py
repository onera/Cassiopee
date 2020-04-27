# - adapts a cells with respect to b points (PyTree) -
# 
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.PyTree as C
import Converter.Internal as I
import Intersector.PyTree as XOR
import KCore.test as test

mesh = G.cart((0,0,0), (1,1,1), (20,20,20))
source = G.cartHexa((8,8,8), (0.2,0.2,0.2), (20,20,20))
#C.convertPyTree2File(mesh, 'm.cgns')
#C.convertPyTree2File(source, 's.cgns')

t = C.newPyTree(['Base',mesh])
t = T.splitNParts(t, 2, multigrid=0, dirs=[1,2,3])

zones = I.getZones(t)
p1 = zones[0]
p1 = C.convertArray2Tetra(p1, split='withBarycenters')
p1 = C.convertArray2NGon(p1)

p2 = C.convertArray2NGon(zones[1])

mesh = XOR.booleanUnion(p1,p2) #conformize the join
#C.convertPyTree2File(mesh, 'u.cgns')

m0 = XOR.adaptCells(mesh,source, sensor_type=0)
m0 = XOR.closeCells(m0)
#C.convertPyTree2File(m0, 'out.cgns')
test.testT(m0,1)

m1 = XOR.adaptCells(mesh,source, sensor_type=0, smoothing_type=1)
m1 = XOR.closeCells(m1)
#C.convertPyTree2File(m1, 'out.cgns')
test.testT(m1,2)

m2 = XOR.adaptCells(mesh,source, sensor_type=1)
m2 = XOR.closeCells(m2)
#C.convertPyTree2File(m2, 'out2.cgns')
test.testT(m2,3)

## dynamic adaptation
hmsh = XOR.createHMesh(mesh)
m = XOR.adaptCells(mesh, source, hmesh = hmsh, sensor_type=0)
m = XOR.conformizeHMesh(m, hmsh)
m = XOR.closeCells(m)
XOR.deleteHMesh(hmsh);
test.testT(m,4)

hmsh = XOR.createHMesh(mesh)
m = XOR.adaptCells(mesh, source, hmesh = hmsh, sensor_type=0, smoothing_type=1)

cm = XOR.conformizeHMesh(m, hmsh)
cm = XOR.closeCells(cm)
test.testT(m,5)

m = XOR.adaptCells(m, source, hmesh = hmsh, sensor_type=0) # applied to existing hmesh with the basic sensor

cm = XOR.conformizeHMesh(m, hmsh)
cm = XOR.closeCells(cm)

XOR.deleteHMesh(hmsh);
test.testT(m,6)
