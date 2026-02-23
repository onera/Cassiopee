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

mesh = C.initVars(mesh, '{centers:Density} = {centers:CoordinateX} + {centers:CoordinateY}')

XOR._setZonesAndJoinsUId(mesh)

m0 = XOR.adaptCells(mesh,source, sensor_type=0)
m0 = XOR.closeCells(m0)

I._rmNodesByName(m0, 'zid')
I._rmNodesByName(m0, 'rid')

test.testT(m0,1)
#C.convertPyTree2File(m0, 'PT_t3_1.cgns')

m1 = XOR.adaptCells(mesh,source, sensor_type=0, smoothing_type=1)
m1 = XOR.closeCells(m1)

I._rmNodesByName(m1, 'zid')
I._rmNodesByName(m1, 'rid')

test.testT(m1,2)
#C.convertPyTree2File(m1, 'PT_t3_2.cgns')

m2 = XOR.adaptCells(mesh,source, sensor_type=0)
m2 = XOR.closeCells(m2)

I._rmNodesByName(m2, 'zid')
I._rmNodesByName(m2, 'rid')

test.testT(m2,3)
#C.convertPyTree2File(m2, 'PT_t3_3.cgns')

## dynamic adaptation
hmsh = XOR.createHMesh(mesh)
m3 = XOR.adaptCells(mesh, source, hmesh=hmsh, sensor_type=0)
m3 = XOR.conformizeHMesh(m3, hmsh)
m3 = XOR.closeCells(m3)

I._rmNodesByName(m3, 'zid')
I._rmNodesByName(m3, 'rid')

XOR.deleteHMesh(hmsh);
test.testT(m3,4)
#C.convertPyTree2File(m3, 'PT_t3_4.cgns')

hmsh = XOR.createHMesh(mesh)
m4 = XOR.adaptCells(mesh, source, hmesh=hmsh, sensor_type=0, smoothing_type=1)

m4 = XOR.conformizeHMesh(m4, hmsh)
m4 = XOR.closeCells(m4)

I._rmNodesByName(m4, 'zid')
I._rmNodesByName(m4, 'rid')

test.testT(m4,5)
#C.convertPyTree2File(m4, 'PT_t3_5.cgns')

m5 = XOR.adaptCells(m4, source, hmesh=hmsh, sensor_type=0) # applied to existing hmesh with the basic sensor

m5 = XOR.conformizeHMesh(m4, hmsh)
m5 = XOR.closeCells(m5)

I._rmNodesByName(m5, 'zid')
I._rmNodesByName(m5, 'rid')

XOR.deleteHMesh(hmsh);
test.testT(m5,6)
#C.convertPyTree2File(m5, 'PT_t3_6.cgns')
