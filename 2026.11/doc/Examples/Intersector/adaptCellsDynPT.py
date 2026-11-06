# - dynamic adaptation
#
import Generator.PyTree as G
import Converter.PyTree as C
import Intersector.PyTree as XOR


mesh = G.cart((0,0,0), (1,1,1), (20,20,20))
mesh = C.convertArray2NGon(mesh)

source = G.cartHexa((8,8,8), (0.2,0.2,0.2), (20,20,20))

XOR._setZonesAndJoinsUId(mesh)

hmsh = XOR.createHMesh(mesh)

senso = XOR.createSensor(hmsh)
XOR.assignData2Sensor(senso, source)

m = XOR.adaptCells(mesh, hmesh=hmsh, sensor=senso)

m = XOR.conformizeHMesh(m, hmsh)
m = XOR.closeCells(m)


XOR.deleteHMesh(hmsh);
XOR.deleteSensor(senso);

C.convertPyTree2File(m, 'out.cgns')
