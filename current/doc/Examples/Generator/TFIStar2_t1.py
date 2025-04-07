# - TFIStar2 (array)
import Converter as C
import Generator as G
import Geom as D
import KCore.test as test

P0 = (0,0,0); P1 = (5,0,0); P2 = (7,3,0); P3 = (4,5,0); P4 = (-2,2,0)

# 5 curves (dont need to be lines)
d1 = D.line(P0, P1, N=11)
d2 = D.line(P1, P2, N=11)
d3 = D.line(P2, P3, N=11)
d4 = D.line(P3, P4, N=11)
d5 = D.line(P4, P0, N=11)

r = G.TFIStar2([d1,d2,d3,d4,d5])
test.testA(r, 1)
