# - close (array) -
import Converter as C
import Generator as G

a = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0.01, 10., (20,20,10))
a = C.convertArray2Tetra(a)
a = G.close(a, 1.e-3)
C.convertArrays2File(a, 'out.plt')
