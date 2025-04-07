# - TFIMono (array) -
import Converter as C
import Generator as G
import Geom as D

a1 = D.circle((0,0,0), 1., tetas=0, tetae=180., N=41)
a2 = D.line((-1,0,0),(1,0,0), N=21)
r = G.TFIMono(a1, a2)
C.convertArrays2File(r, 'out.plt')
