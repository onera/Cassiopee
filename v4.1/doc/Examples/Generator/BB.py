# - BB (array) -
import Generator as G
import Converter as C
import Geom as D

s = D.circle((0,0,0), 1., N=100)
a = G.BB(s)
C.convertArrays2File([a,s], 'out.plt')
