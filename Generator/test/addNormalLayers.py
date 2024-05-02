# - addNormalLayers (array) -
import Generator as G
import Converter as C
import Geom as D

d = C.array('d', 3, 1, 1)
d[1][0,0] = 0.1; d[1][0,1] = 0.2; d[1][0,2] = 0.3
a = D.sphere( (0,0,0), 1, 50 )
a = G.addNormalLayers(a, d, niter=20)
C.convertArrays2File([a], 'out.plt')
