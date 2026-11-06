# - splitBAR (array) -
import Generator as G
import Converter as C
import Transform as T

a = G.cart((0,0,0), (1,1,1), (50,1,1))
a = C.convertArray2Tetra(a)
b = T.splitBAR(a, 5)
C.convertArrays2File(b, 'out.plt')
