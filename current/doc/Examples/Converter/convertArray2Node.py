# - convertArray2Node (array) -
import Converter as C
import Generator as G

a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
b = C.convertArray2Node(a)
C.convertArrays2File([b], 'out.plt')
