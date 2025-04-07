# - TTM (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Geom.PyTree as D

P0 = (0,0,0); P1 = (5,0,0); P2 = (0,7,0); P3 = (5,7,0)

# Geometry
d1 = D.line(P0,P1); d2 = D.line(P2,P3)
d3 = D.line(P0,P2); d4 = D.line(P1,P3)

# Regular discretisation of each line
Ni = 20; Nj = 10
r = G.cart((0,0,0), (1./(Ni-1),1,1), (Ni,1,1))
q = G.cart((0,0,0), (1./(Nj-1),1,1), (Nj,1,1))
r1 = G.map(d1, r); r2 = G.map(d2, r)
r3 = G.map(d3, q); r4 = G.map(d4, q)

# TTM
m = G.TFI([r1, r2, r3, r4])
m = G.TTM(m)
t = C.newPyTree(['Base',2]); t[2][1][2].append(m)
C.convertPyTree2File(t, 'out.cgns')
