# - densify (array) -
import Generator as G
import Converter as C
import Geom as D

a = D.circle((0,0,0), 1., 10)
b = G.densify(a, 0.01)
C.convertArrays2File(b, 'out.plt')
