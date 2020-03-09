# - adapts a cells with respect to b points (PyTree) -
import Intersector.PyTree as XOR
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.1), (5,5,5))
a = C.convertArray2NGon(a); a = G.close(a)
#C.convertPyTree2File([a], 'a.cgns')
b = G.cartHexa((0.,0.,0.), (0.005,0.005,0.005), (5,5,5))
#C.convertPyTree2File([b], 'b.cgns')

a = C.fillEmptyBCWith(a, 'wall', 'BCWall')

m = XOR.adaptCells(a,b, sensor_type=0)
m = XOR.closeCells(m)
C.convertPyTree2File(m, 'out.cgns')

m = XOR.adaptCells(a,b, sensor_type=1)
m = XOR.closeCells(m)
C.convertPyTree2File(m, 'xout.cgns')

