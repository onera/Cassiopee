# - splitSize (array) -
import Generator as G
import Transform as T
import Converter as C

a = G.cart((0,0,0),(1,1,1),(50,20,10))
B = T.splitSize(a, 2000, type=0)
C.convertArrays2File(B, 'out.plt')
