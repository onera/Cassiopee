# - lineDrive (pyTree)-
import Geom.PyTree as D
import Converter.PyTree as C

# With one driving curve
a = D.naca(12.)
l = D.line((0,0,0), (0,0.,1.))
a = D.lineDrive(a, l)
C.convertPyTree2File(a, 'out.cgns')

# With a set of driving curves
a = D.naca(12.)
d1 = D.line((0,0,0), (0.,0.,1.))
d2 = D.line((1,0,0), (2,0,1))
a = D.lineDrive(a, [d1,d2])
C.convertPyTree2File(a, 'out.cgns')
