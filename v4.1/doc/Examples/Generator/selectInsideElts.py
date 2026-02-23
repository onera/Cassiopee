# - selectInsideElts (array) -
import Converter as C
import Generator as G
import Geom as D

a = G.cart( (0,0,0), (1,1,1), (10,10,1)); a = C.convertArray2Tetra(a)
b = D.circle( (5,5,0), 3.); b = C.convertArray2Tetra(b)
a = G.selectInsideElts(a, [b])
C.convertArrays2File([a,b], 'out.plt')
