# - getCircumCircleMap (pyTree) -
import Geom.PyTree as D
import Generator.PyTree as G
import Converter.PyTree as C

a = D.sphere((0,0,0), 1, 50)
a = C.convertArray2Tetra(a)
t = C.newPyTree(['Base',2,a])
t = G.getCircumCircleMap(t)
C.convertPyTree2File(t, 'out.cgns')
