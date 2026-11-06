# - lineDrive (array) -
import Geom as D
import Converter as C

# With one driving curve
a = D.naca(12.)
b = D.line((0,0,0), (0.,0.,1.))
c = D.lineDrive(a, b)
C.convertArrays2File([c], 'out.plt')

# With a set of driving curves
a = D.naca(12.)
d1 = D.line((0,0,0), (0.,0.,1.))
d2 = D.line((1,0,0), (2,0,1))
c = D.lineDrive(a, [d1,d2])
C.convertArrays2File([c,d1,d2,a], 'out.plt')
