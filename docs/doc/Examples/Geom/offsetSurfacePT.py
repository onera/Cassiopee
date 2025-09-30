# - offsetSurface (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C
import Geom.Offset

a = D.circle((0,0,0), 1.)
b = Geom.Offset.offsetSurface(a, offset=1., pointsPerUnitLength=10., algo=0, dim=2)
C.convertPyTree2File(a, 'out.cgns')
