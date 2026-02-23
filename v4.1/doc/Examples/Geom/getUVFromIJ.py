# - getUVFromIJ (array) -
import Geom as D
import Generator as G
import Converter as C

a = G.cart((0,0,0), (1,1,1), (10,10,1))

a = D.getUVFromIJ(a)

C.convertArrays2File(a, 'out.plt')
