# - adapts a cells with respect to b points (PyTree) -
import Intersector.PyTree as XOR
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test
import Converter.Internal as I

a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.1), (5,5,5))
a = C.convertArray2NGon(a); a = G.close(a)
#C.convertPyTree2File(a, 'a.cgns')
b = G.cartHexa((0.,0.,0.), (0.005,0.005,0.005), (5,5,5))
#C.convertPyTree2File(b, 'b.cgns')

a = C.fillEmptyBCWith(a, 'wall', 'BCWall')
a = C.initVars(a, '{centers:Density} = {centers:CoordinateX} + {centers:CoordinateY}')
#C.convertPyTree2File(a, 'm0.cgns')

XOR._setZonesAndJoinsUId(a)

# static adaptation
m = XOR.adaptCells(a,b, sensor_type=0)
m = XOR.closeCells(m)

I._rmNodesByName(m, 'zid')
I._rmNodesByName(m, 'rid')

test.testT(m,1)
#C.convertPyTree2File(m, 'PT_t1_1.cgns')

m = XOR.adaptCells(a,b, sensor_type=1, smoothing_type=1)
m = XOR.closeCells(m)

I._rmNodesByName(m, 'zid')
I._rmNodesByName(m, 'rid')

test.testT(m,2)
#C.convertPyTree2File(m, 'PT_t1_2.cgns')

m = XOR.adaptCells(a,b, sensor_type=0, smoothing_type=1)
m = XOR.closeCells(m)

I._rmNodesByName(m, 'zid')
I._rmNodesByName(m, 'rid')

test.testT(m,3)
#C.convertPyTree2File(m, 'PT_t1_3.cgns')

# dynamic adaptation
hmsh = XOR.createHMesh(a)
m = XOR.adaptCells(a, b, hmesh=hmsh, sensor_type=0)
m = XOR.conformizeHMesh(m, hmsh)
m = XOR.closeCells(m)

I._rmNodesByName(m, 'zid')
I._rmNodesByName(m, 'rid')

XOR.deleteHMesh(hmsh);
test.testT(m,4)
#C.convertPyTree2File(m, 'PT_t1_4.cgns')

hmsh = XOR.createHMesh(a)
m = XOR.adaptCells(a, b, hmesh=hmsh, sensor_type=0, smoothing_type=1)

cm = XOR.conformizeHMesh(m, hmsh)
cm = XOR.closeCells(cm)

I._rmNodesByName(cm, 'zid')
I._rmNodesByName(cm, 'rid')

test.testT(cm,5)
#C.convertPyTree2File(cm, 'PT_t1_5.cgns')

m = XOR.adaptCells(m, b, hmesh=hmsh, sensor_type=0) # applied to existing hmesh with the basic sensor

cm = XOR.conformizeHMesh(cm, hmsh)
cm = XOR.closeCells(cm)

I._rmNodesByName(cm, 'zid')
I._rmNodesByName(cm, 'rid')

XOR.deleteHMesh(hmsh);
test.testT(cm,6)
#C.convertPyTree2File(cm, 'PT_t1_6.cgns')
