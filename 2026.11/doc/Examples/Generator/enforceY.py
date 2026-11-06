# - enforceY (array) -
import Generator as G
import Converter as C

# Distribution
Ni = 50; Nj = 50
a = G.cart((0,0,0), (1./(Ni-1), 0.5/(Nj-1),1), (Ni,Nj,1))
b = G.enforceY(a, 0.3, 0.001, (10,15))
C.convertArrays2File([b], "out.plt")
