# - enforceMoinsZ (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C

Ni = 50; Nj = 2; Nk = 50
a = G.cart((0,0,0), (1./(Ni-1), 1., 0.5/(Nk-1)), (Ni,Nj,Nk))
b = G.enforceMoinsZ(a, 1.e-3, 10,15)
C.convertPyTree2File(b, 'out.cgns')
