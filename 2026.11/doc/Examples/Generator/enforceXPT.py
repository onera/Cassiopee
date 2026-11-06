# - enforceX (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C

Ni = 50; Nj = 50; Nk = 2
a = G.cart((0,0,0), (1./(Ni-1), 0.5/(Nj-1),1), (Ni,Nj,Nk))
a = G.enforceX(a, 0.3, 0.001, (13,25))
C.convertPyTree2File(a, 'out.cgns')
