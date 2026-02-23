# - enforcePlusZ (array) -
import Generator as G
import Converter as C

Ni = 50; Nj = 1; Nk = 50
a = G.cart((0,0,0), (1./(Ni-1), 1., 0.5/(Nk-1)), (Ni,Nj,Nk))
b = G.enforcePlusZ(a, 1.e-3, 10,20)
C.convertArrays2File([b], "out.plt")
