# - adapts a cells with respect to b points (PyTree) -
import Intersector.PyTree as XOR
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.1), (5,5,5))
a = C.convertArray2NGon(a); a = G.close(a)
#C.convertPyTree2File(a, 'a.cgns')
b = G.cartHexa((0.,0.,0.), (0.005,0.005,0.005), (5,5,5))
#C.convertPyTree2File(b, 'b.cgns')

a = C.fillEmptyBCWith(a, 'wall', 'BCWall')
a = C.initVars(a, '{centers:Density} = {centers:CoordinateX} + {centers:CoordinateY}')

XOR._setZonesAndJoinsUId(a)

## static adaptation
m = XOR.adaptCells(a,b, sensor_type=0)
m = XOR.closeCells(m)
C.convertPyTree2File(m, 'out.cgns')

m = XOR.adaptCells(a,b, sensor_type=1)
m = XOR.closeCells(m)
C.convertPyTree2File(m, 'xout.cgns')

m = XOR.adaptCells(a,b, sensor_type=0, smoothing_type=1)
m = XOR.closeCells(m)
C.convertPyTree2File(m, 'out2.cgns')

## dynamic adaptation
hmsh = XOR.createHMesh(a)
m = XOR.adaptCells(a, b, hmesh=hmsh, sensor_type=0)
cm = XOR.conformizeHMesh(m, hmsh)
cm = XOR.closeCells(m)
XOR.deleteHMesh(hmsh);
C.convertPyTree2File(cm, 'out3.cgns')

hmsh = XOR.createHMesh(a)
m = XOR.adaptCells(a, b, hmesh=hmsh, sensor_type=0, smoothing_type=1)

cm = XOR.conformizeHMesh(m, hmsh)
cm = XOR.closeCells(cm)
C.convertPyTree2File(cm, 'out4.cgns')

m = XOR.adaptCells(m, b, hmesh=hmsh, sensor_type=0) # applied to existing hmesh with the geometrical sensor

cm = XOR.conformizeHMesh(cm, hmsh)
cm = XOR.closeCells(cm)

XOR.deleteHMesh(hmsh);
C.convertPyTree2File(cm, 'out5.cgns')
