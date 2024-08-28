# - splitCurvatureAngle (array) -
import Converter as C
import Transform as T
import Geom as D

line = D.line((0.,0.,0.), (1.,1.,0.), 51)
line2 = D.line((-1.,1.,0.), (0.,0.,0.), 51)
a = T.join(line2, line)
list = T.splitCurvatureAngle(a, 30.)
C.convertArrays2File(list, 'out.plt')
