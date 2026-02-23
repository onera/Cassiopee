# - TFITri (array)
import Converter as C
import Generator as G
import Geom as D

P0 = (0,0,0); P1 = (5,0,0); P2 = (1,7,0)

# 3 curves (dont need to be lines)
d1 = D.line(P0, P1, N=11)
d2 = D.line(P1, P2, N=11)
d3 = D.line(P0, P2, N=11)
r = G.TFITri(d1, d2, d3)
C.convertArrays2File(r, 'out.plt')
