# - convertArray2Hexa (array) -
import Converter as C
import Generator as G

# 2D: quad
a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
a = C.convertArray2Hexa(a)

# 3D: hexa
b = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
b = C.convertArray2Hexa(b)
C.convertArrays2File([a,b], 'out.plt')
