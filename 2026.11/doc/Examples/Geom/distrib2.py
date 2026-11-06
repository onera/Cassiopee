# - distrib2 (array) -
import Geom as D
import Generator as G
import Converter as C

a = D.line((0,0,0), (4,4,0), N=30)
b = D.distrib2(a, 0.001, 0.1)
c = G.map(a, b)
C.convertArrays2File(c, 'out.plt')
