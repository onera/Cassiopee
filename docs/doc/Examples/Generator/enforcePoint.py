# - enforcePoint (array) -
import Converter as C
import Generator as G

# distribution
Ni = 20; Nj = 20
a = G.cart((0,0,0), (1./(Ni-1),5./(Nj-1),1), (Ni,Nj,1))
b = G.enforcePoint(a, 0.5)
C.convertArrays2File([b], "out.plt")
