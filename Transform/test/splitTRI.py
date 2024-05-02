# - splitTRI (array) -
import Generator as G
import Converter as C
import Transform as T

a = G.cart((0,0,0),(1,1,1),(5,5,1))
a = C.convertArray2Tetra(a)
C.convertArrays2File(a, 'out.plt')
c = [[10,16,22], [2,8,9]]
d = T.splitTRI(a, c)

C.convertArrays2File(d[0], 'out1.plt')
C.convertArrays2File(d[1], 'out2.plt')
