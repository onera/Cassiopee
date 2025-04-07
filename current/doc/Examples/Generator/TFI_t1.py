# - TFI 2D structure -
import Generator as G
import Geom as D
import Converter as C
import KCore.test as test

P0 = (0,0,0)
P1 = (5,0,0)
P2 = (0,7,0)
P3 = (5,7,0)

# Geometrie
d1 = D.line(P0, P1)
d2 = D.line(P2, P3)
pts = C.array('x,y,z', 5, 1, 1)
x = pts[1][0]; y = pts[1][1]; z = pts[1][2]

x[0 ] = 0. ; y[0] = 0.;
x[1 ] =-2. ; y[1] = 2.;
x[2 ] =-3. ; y[2] = 3.;
x[3 ] = 2. ; y[3] = 5.;
x[4 ] = 0. ; y[4] = 7.;
b1 = D.bezier(pts)
#
pts = C.array('x,y,z', 5, 1, 1)
x = pts[1][0]; y = pts[1][1]; z = pts[1][2]
x[0 ] = 5.; y[ 0 ] = 0.;
x[1 ] = 3.; y[ 1 ] = 2.;
x[2 ] = 2.; y[ 2 ] = 3.;
x[3 ] = 6.; y[ 3 ] = 5.;
x[4 ] = 5.; y[ 4 ] = 7.;

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
test.testA([m],1)
