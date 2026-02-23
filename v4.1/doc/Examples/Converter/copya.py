# - copy (array) -
import Converter as C
import Generator as G

a = G.cart((0,0,0,), (1,1,1), (10,10,10))
b = C.copy(a)
C.convertArrays2File([b], "out.plt")
