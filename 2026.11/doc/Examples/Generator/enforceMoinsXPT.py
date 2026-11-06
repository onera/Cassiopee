# - enforceMoinsX (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

Ni = 50; Nj = 50; Nk = 1
a = G.cart((0,0,0), (1./(Ni-1), 0.5/(Nj-1),1), (Ni,Nj,Nk))
b = G.enforceMoinsX(a, 1.e-3, (10,15))
C.convertPyTree2File(b, 'out.cgns')
