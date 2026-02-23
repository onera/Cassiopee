# - convertArray2NGon(array) -
import Converter as C
import Generator as G

a = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (2,2,1))
b = C.convertArray2NGon(a)
C.convertArrays2File(b, 'out.plt')
