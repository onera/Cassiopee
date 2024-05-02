# - text1D (array) -
import Geom as D
import Converter as C
import Transform as T

a = D.text1D("Cassiopee - text1")
b = D.text1D("Cassiopee - text1 smoothed", smooth=4, offset=1.)
b = T.translate(b, (0,-12,0))
c = D.text1D("Cassiopee - vera", font='vera')
c = T.translate(c, (0,-24,0))
d = D.text1D("Cassiopee - chancery", font='chancery')
d = T.translate(d, (0,-36,0))
e = D.text1D("Cassiopee - courier", font='courier')
e = T.translate(e, (0,-48,0))
f = D.text1D("Cassiopee - nimbus", font='nimbus')
f = T.translate(f, (0,-60,0))

C.convertArrays2File(a+b+c+d+e+f, 'out.plt')
