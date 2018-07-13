# - enforceLine (pyTree)-
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D

Ni = 50; Nj = 50
a = G.cart((0,0,0), (1./(Ni-1), 0.5/(Nj-1),1), (Ni,Nj,1))
b = D.line((0.,0.2,0.), (1.,0.2,0.), 20)
a = G.enforceLine(a, b, 0.01, (5,3))
C.convertPyTree2File(a, 'out.cgns')
