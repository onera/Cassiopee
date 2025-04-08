# - smooth (array) -
import Transform as T
import Converter as C
import Geom as D
a = D.sphere6((0,0,0), 1, N=20)
b = T.smooth(a, eps=0.5, niter=20)
C.convertArrays2File(a+b, "out.plt")
