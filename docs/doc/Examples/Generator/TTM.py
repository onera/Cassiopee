# - TTM (array) -
import Converter as C
import Generator as G
import Geom as D

P0 = (0,0,0); P1 = (5,0,0); P2 = (0,7,0); P3 = (5,7,0)

# Geometry
d1 = D.line(P0, P1); d2 = D.line(P2, P3)
pts = C.array('x,y,z', 5, 1, 1)
x = pts[1][0]; y = pts[1][1]; z = pts[1][2]

x[0] = 0. ; y[0] = 0.; z[0] = 0.
x[1] =-2. ; y[ 1 ] = 2.; z[1] = 0.
x[2] =-3. ; y[ 2 ] = 3.; z[2] = 0.
x[3] = 2. ; y[ 3 ] = 5.; z[3] = 0.
x[4] = 0. ; y[ 4 ] = 7.; z[4] = 0.
b1 = D.bezier(pts)

x[0] = 5.; y[ 0 ] = 0.; z[ 0 ] = 0.
x[1] = 3.; y[ 1 ] = 2.; z[ 1 ] = 0.
x[2] = 2.; y[ 2 ] = 3.; z[ 2 ] = 0.
x[3] = 6.; y[ 3 ] = 5.; z[ 3 ] = 0.
x[4] = 5.; y[ 4 ] = 7.; z[ 4 ] = 0.
b2 = D.bezier( pts )
C.convertArrays2File([d1, d2, b1, b2], 'geom.plt')

# Regular discretision of each line
Ni = 20; Nj = 10
r = G.cart((0,0,0), (1./(Ni-1),1,1), (Ni,1,1))
q = G.cart((0,0,0), (1./(Nj-1),1,1), (Nj,1,1))
r1 = G.map(d1, r)
r2 = G.map(d2, r)
r3 = G.map(b1, q)
r4 = G.map(b2, q)

# TTM
m = G.TFI([r1, r2, r3, r4])
m2 = G.TTM(m, 2000)
C.convertArrays2File([m,m2], 'out.plt')
