# - enforcePoint (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

Ni = 20; Nj = 20
a = G.cart((0,0,0), (1./(Ni-1),5./(Nj-1),1), (Ni,Nj,1))
b = G.enforcePoint(a, 0.5)
C.convertPyTree2File(b, 'out.cgns')
