# - convertHexa2Struct (array) -
import Converter as C
import Generator as G

# 2D:
a = G.cart((0.,0.,0.), (1,1,2), (5,5,1))
a = C.convertArray2Hexa(a)
b = C.convertHexa2Struct(a)
print(b)
C.convertArrays2File([b], "out.plt")

# 3D:
