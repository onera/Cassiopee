# - TFIO (array) -
import Converter as C
import Generator as G
import Geom as D

a = D.circle((0,0,0), 1., N=41)
r = G.TFIO(a)
C.convertArrays2File(r, 'out.plt')
