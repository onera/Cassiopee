# - TFI (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Geom.PyTree as D

# Geometry
P0 = (0,0,0); P1 = (5,0,0); P2 = (0,7,0); P3 = (5,7,0)
Ni = 20; Nj = 10
d1 = D.line(P0, P1,Ni); d2 = D.line(P2, P3,Ni)
d3 = D.line(P0, P2,Nj); d4 = D.line(P1, P3,Nj)
m = G.TFI([d1, d2, d3, d4])
C.convertPyTree2File(m, 'out.cgns')
