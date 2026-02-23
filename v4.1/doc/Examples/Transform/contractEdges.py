# - contractEdges (array) -
import Geom as D
import Converter as C
import Transform as T
import Generator as G


a = D.triangle((0,0,0), (1,0,0), (1,1,0))
b = D.triangle((1,0,0), (2,0,0), (1,0.5,0))
c = D.triangle((1,0,0), (1,1,0), (1,0.5,0))

r = T.join([a,b,c])
r = G.close(r)

r = T.transform.contractEdges(r, 1)

C.convertArrays2File(r, 'out.plt')
