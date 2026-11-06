# - convertTri2Quad (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
a, b = C.convertTri2Quad(a, 30.)
C.convertPyTree2File([a,b], 'out.cgns')
