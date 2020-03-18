# - adapts a cells with respect to b points (PyTree) -
import Intersector.PyTree as XOR
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.1), (5,5,5))
a = C.convertArray2Tetra(a, split='withBarycenters')
a = C.convertArray2NGon(a); a = G.close(a)
#C.convertPyTree2File(a, 'a.cgns')
b = G.cartHexa((0.,0.,0.), (0.005,0.005,0.005), (5,5,5))
#C.convertPyTree2File(b, 'b.cgns')

a = C.fillEmptyBCWith(a, 'wall', 'BCWall')

## static adaptation
m = XOR.adaptCells(a,b, sensor_type=0)
m = XOR.closeCells(m)
test.testT(m,1)

m = XOR.adaptCells(a,b, sensor_type=1)
m = XOR.closeCells(m)
test.testT(m,2)

m = XOR.adaptCells(a,b, sensor_type=2)
m = XOR.closeCells(m)
test.testT(m,3)

## dynamic adaptation
hmsh = XOR.createHMesh(a)
m = XOR.adaptCells(a, b, hmesh = hmsh, sensor_type=0)
m = XOR.conformizeHMesh(m, hmsh)
m = XOR.closeCells(m)
XOR.deleteHMesh(hmsh);
test.testT(m,4)

hmsh = XOR.createHMesh(a)
m = XOR.adaptCells(a, b, hmesh = hmsh, sensor_type=2)

cm = XOR.conformizeHMesh(m, hmsh)
cm = XOR.closeCells(cm)
test.testT(m,5)

m = XOR.adaptCells(m, b, hmesh = hmsh, sensor_type=0) # applied to existing hmesh with the basic sensor

cm = XOR.conformizeHMesh(m, hmsh)
cm = XOR.closeCells(cm)

XOR.deleteHMesh(hmsh);
test.testT(m,6)

