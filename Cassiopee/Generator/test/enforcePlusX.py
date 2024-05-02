# - enforcePlusX (array) -
import Generator as G
import Converter as C

# Distribution
Ni = 50; Nj = 50
a = G.cart((0,0,0), (1./(Ni-1), 0.5/(Nj-1),1), (Ni,Nj,1))
b = G.enforcePlusX(a, 1.e-3, (10,20))
C.convertArrays2File([b], "out.plt")
