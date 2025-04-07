# - TFI 2D structured (array)
# - TFI 3D structured (array)
# - TFI TRI (array)
import Converter as C
import Generator as G
import Geom as D
import Transform as T

#------------------
# TFI 2D structured
#------------------
P0 = (0,0,0); P1 = (5,0,0); P2 = (0,7,0); P3 = (5,7,0)

# Geometrie
d1 = D.line(P0, P1); d2 = D.line(P2, P3)

pts = C.array('x,y,z',5,1,1)
x = pts[1][0]; y = pts[1][1]; z = pts[1][2]
x[0] = 0.; y[0] = 0.; z[0] = 0.
x[1] =-2.; y[1] = 2.; z[1] = 0.
x[2] =-3.; y[2] = 3.; z[2] = 0.
x[3] = 2.; y[3] = 5.; z[3] = 0.
x[4] = 0.; y[4] = 7.; z[4] = 0.
b1 = D.bezier(pts)

pts = C.array('x,y,z',5,1,1)
x = pts[1][0]; y = pts[1][1]; z = pts[1][2]
x[0] = 5.; y[0] = 0.; z[0] = 0.
x[1] = 3.; y[1] = 2.; z[1] = 0.
x[2] = 2.; y[2] = 3.; z[2] = 0.
x[3] = 6.; y[3] = 5.; z[3] = 0.
x[4] = 5.; y[4] = 7.; z[4] = 0.
b2 = D.bezier( pts )

C.convertArrays2File([d1, d2, b1, b2], "geom.plt")

# Regular discretisation of each line
Ni = 20
Nj = 10
r = G.cart((0,0,0), (1./(Ni-1),1,1), (Ni,1,1))
q = G.cart((0,0,0), (1./(Nj-1),1,1), (Nj,1,1))
r1 = G.map(d1, r)
r2 = G.map(d2, r)
r3 = G.map(b1, q)
r4 = G.map(b2, q)

# TFI 2D
m = G.TFI([r1, r2, r3, r4])
C.convertArrays2File([r1,r2,r3,r4,m], 'tfi2d.plt')

#------------------
# TFI 3D structured
#------------------
xo = 0.; yo = 0.; zo = 0.
nx = 21; ny = 21; nz = 21
hx = 1./(nx-1); hy = 1./(ny-1); hz = 1./(nz-1)

# z = cste
fzmin = G.cart((xo,yo,zo), (hx,hy,1.), (nx,ny,1))
fzmax = T.translate(fzmin, (0.,0.,1.))

# x = cste
fxmin = G.cart((xo,yo,zo),(1,hy,hz),(1,ny,nz))
fxmin = T.reorder(fxmin,(3,1,2))
fxmax = T.translate(fxmin, (1.,0.,0.))

# y = cste
fymin = G.cart((xo,yo,zo),(hx,1.,hz),(nx,1,nz))
fymin = T.reorder(fymin,(1,3,2))
fymax = T.translate(fymin, (0.,1.,0.))

r = [fxmin,fxmax,fymin,fymax,fzmin,fzmax]
m = G.TFI(r)
C.convertArrays2File(r+[m], 'tfi3d.plt')

#---------
# TFI TRI
#---------
l1 = D.line((0,0,0),(0,1,0), 15)
l2 = D.line((0,0,0),(1,0,0), 15)
l3 = D.line((1,0,0),(0,1,0), 15)
tri = G.TFI([l1,l2,l3])
C.convertArrays2File([tri], 'tfitri.plt')
