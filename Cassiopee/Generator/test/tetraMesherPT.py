# - tetraMesher (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Post.PyTree as P
import Transform.PyTree as T

a = G.cart((0,0,0), (1,1,1), (3,3,3))
ext = P.exteriorFaces(a)
ext = C.convertArray2Tetra(ext)
ext = G.close(ext)
ext = T.reorder(ext, (-1,))
m = G.tetraMesher(ext, algo=1)
C.convertPyTree2File(m, 'out.cgns')
