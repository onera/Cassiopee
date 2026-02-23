# - gapsmanager (array) -
import Geom as D
import Converter as C
import Generator as G

a = D.sphere6((0,0,0), 1, N=10)
a = C.node2Center(a)
a = C.convertArray2Tetra(a)
b = G.gapsmanager(a, mode=2)
C.convertArrays2File(b, 'out.plt')
