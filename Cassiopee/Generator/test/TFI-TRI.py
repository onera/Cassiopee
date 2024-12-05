# - TFI TRI (array) -
import Generator as G
import Converter as C
import Geom as D

imin = D.line((0,0,0), (0,1,0), 15)
jmin = D.line((0,0,0), (1,0,0), 15)
diag = D.line((1,0,0), (0,1,0), 15)

pts = C.array('x,y,z', 5, 1, 1)
x = pts[1][0]; y = pts[1][1]; z = pts[1][2]

x[0] = 0.;  y[0] = 0.;
x[1] = 0.2; y[1] =-0.2;
x[2] = 0.4; y[2] =-0.5;
x[3] = 0.8; y[3] =-0.1;
x[4] = 1.;  y[4] = 0.;
b2 = D.bezier(pts, 15)

N = 5
c = G.cart((0,0,0),(1./(N-1),1,1),(N,1,1))
c = G.enforceX(c, 0.5, 0.01, (2, 5))
b2 = G.map(b2, c)

C.convertArrays2File([imin,b2,diag], 'contours.plt', 'bin_tp')

tri = G.TFI([imin, b2, diag])
C.convertArrays2File([tri], "out.plt", "bin_tp")
