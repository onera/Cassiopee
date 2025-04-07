# - offsetSurface (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C

a = D.circle((0,0,0), 1.)
b = D.offsetSurface(a, offset=1., pointsPerUnitLength=10., algo=0, dim=2)
C.convertPyTree2File(a, 'out.plt')
