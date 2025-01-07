# - TTM (pyTree) -
import Generator.PyTree as G
import Geom.PyTree as D
import KCore.test as test

P0 = (0,0,1)
P1 = (5,0,1)
P2 = (0,7,1)
P3 = (5,7,1)

# Geometrie
d1 = D.line(P0, P1)
d2 = D.line(P2, P3)

P0 = ( 0,0,1); P1 = (-2,2,1); P2 = (-3,3,1); P3 = (2,5,1); P4 = ( 0,7,1)
pts = D.polyline([P0,P1,P2,P3,P4])
b1 = D.bezier(pts)

P0 = (5,0,1); P1 = (3,2,1); P2 = (2,3,1); P3 = (6,5,1); P4 = (5,7,1)
pts = D.polyline([P0,P1,P2,P3,P4])
b2 = D.bezier(pts)

# Discretisation reguliere de chaque ligne
Ni = 20; Nj = 10
r = G.cart((0,0,0), (1./(Ni-1),1,1), (Ni,1,1))
q = G.cart((0,0,0), (1./(Nj-1),1,1), (Nj,1,1))
r1 = G.map(d1, r)
r2 = G.map(d2, r)
r3 = G.map(b1, q)
r4 = G.map(b2, q)

# TTM
m = G.TFI([r1, r2, r3, r4])
m = G.TTM(m)
test.testT(m,1)
