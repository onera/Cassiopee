# - addSeparationLine (pyTree)-
import Geom.PyTree as D
import Converter.PyTree as C

a1 = D.circle((0,0,0), 1, 0., 360, 1000)
a2 = D.line((0.,1.,0.), (0.,2.,0), 100)
zones = D.addSeparationLine(a1, a2)
C.convertPyTree2File(zones, 'out.cgns')
