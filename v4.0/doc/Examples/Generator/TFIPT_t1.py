# - TFI 2D structure -
import Generator.PyTree as G
import Geom.PyTree as D
import KCore.test as test

P0 = (0,0,0); P1 = (5,0,0); P2 = (0,7,0); P3 = (5,7,0)

# Geometrie
d1 = D.line(P0, P1); d2 = D.line(P2, P3)
pts = D.polyline([(0.,0.,0.),(-2,2,0),(-3,3,0.),(2,5,0.),(0.,7,0.)])
b1 = D.bezier(pts)
pts = D.polyline([(5.,0.,0.),( 3,2,0),( 2,3,0.),(6,5,0.),(5.,7.,0.)])
b2 = D.bezier( pts )

# Discretisation reguliere de chaque ligne
Ni = 20; Nj = 10
r = G.cart((0,0,0), (1./(Ni-1),1,1), (Ni,1,1))
q = G.cart((0,0,0), (1./(Nj-1),1,1), (Nj,1,1))
r1 = G.map(d1, r)
r2 = G.map(d2, r)
r3 = G.map(b1, q)
r4 = G.map(b2, q)

# TFI
r = [r1,r2,r3,r4]
m = G.TFI(r)
test.testT(m)
