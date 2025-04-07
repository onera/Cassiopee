# - signNGonFaces (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cartNGon((0,0,0), (1,1,1), (3,3,2))
a = C.signNGonFaces(a)
C.convertPyTree2File(a, 'out.cgns')
