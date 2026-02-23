# - plaster (array) -
import Generator as G
import Converter as C
import Transform as T
import Post as P
import Geom as D

a = D.sphere( (0,0,0), 1, N=30)
a = T.subzone(a, (6,1,1), (a[2]//2,a[3],a[4]))
a = C.convertArray2Hexa(a); a = G.close(a)

# contours
c = P.exteriorFaces(a)
cs = T.splitConnexity(c)

# plaster hole
p = G.plaster([cs[0]], [a])
C.convertArrays2File([a,p]+cs, 'out.plt')
