# - convertTri2Quad (array) -
import Converter as C
import Generator as G

a = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
a, b = C.convertTri2Quad(a, 30.)
C.convertArrays2File([a,b], 'out.plt')
