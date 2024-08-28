# - enforcePlusY (array) -
import Generator as G
import Converter as C

Ni = 50; Nj = 50
a = G.cart((0,0,0), (1./(Ni-1), 0.5/(Nj-1),1), (Ni,Nj,1))
b = G.enforcePlusY(a, 1.e-3, (10,15))
C.convertArrays2File([b], "out.plt")
