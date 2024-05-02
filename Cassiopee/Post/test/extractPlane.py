# - extractPlane (array) -
import Converter as C
import Post as P
import Transform as T
import Generator as G

m = G.cylinder((0,0,0), 1, 5, 0., 360., 10., (50,50,50))
m = T.rotate(m , (0,0,0), (1,0,0), 35.)
a = P.extractPlane([m], (0.5, 1., 0., 1), 2)
C.convertArrays2File([m,a], "out.plt")
