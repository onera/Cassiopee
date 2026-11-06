# - tetraMesher (array) -
import Generator as G
import Converter as C
import Post as P
import Transform as T

a = G.cart((0,0,0), (1,1,1), (3,3,3))
ext = P.exteriorFaces(a)
ext = C.convertArray2Tetra(ext)
ext = G.close(ext)
ext = T.reorder(ext, (-1,))
m = G.tetraMesher(ext, algo=1)
C.convertArrays2File(m, 'out.plt')
