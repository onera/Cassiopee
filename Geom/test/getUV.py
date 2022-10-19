# - getUV (array) -
import Geom as D
import Generator as G
import Converter as C
import Post as P

a = G.cartTetra((0,0,0), (1,1,1), (10,10,10))
a = P.exteriorFaces(a)
a = C.initVars(a, 'u', 0.)
a = C.initVars(a, 'v', 0.)
ret = D.getUV(a, 2., 0.)

C.convertArrays2File(ret[0], 'out.plt')
C.convertArrays2File(ret[1], 'color.png')
C.convertArrays2File(ret[2], 'bump.png')
