# - getTriQualitylityMap (PyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D

a = D.sphere((0,0,0), 1, N=10)
a = C.convertArray2Tetra(a)
a = G.close(a)
t = C.newPyTree(['Base',2,a])
t = G.getTriQualityMap(t)
C.convertPyTree2File(t, 'out.cgns')
