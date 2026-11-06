# - enforceY (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C

Ni = 50; Nj = 50
a = G.cart((0,0,0), (1./(Ni-1), 0.5/(Nj-1),1), (Ni,Nj,1))
b = G.enforceY(a, 0.3, 0.001, (10,15))
C.convertPyTree2File(b, 'out.cgns')
