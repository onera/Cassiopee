# - mapCurvature (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D

ni = 2; nj = 3
a = G.cart((0,0,0), (1,1,1), (ni,nj,1))
C.setValue(a,'GridCoordinates', (1,1,1), [1.,1.,2.])
C.setValue(a,'GridCoordinates', (1,2,1), [1.,2.,5.])
C.setValue(a,'GridCoordinates', (1,3,1), [1.,3.,2.])
C.setValue(a,'GridCoordinates', (2,1,1), [2.,1.,2.])
C.setValue(a,'GridCoordinates', (2,2,1), [2.,2.,5.])
C.setValue(a,'GridCoordinates', (2,3,1), [2.,3.,2.])
b = D.bezier(a, density=10.)

b = G.mapCurvature(b, N=100, power=0.5, dir=1)

C.convertPyTree2File(b, 'out.cgns')
