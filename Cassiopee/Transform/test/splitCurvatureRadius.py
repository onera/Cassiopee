# - splitCurvatureRadius (array) -
import Converter as C
import Transform as T
import Geom as D

pts = C.array('x,y,z', 7, 1, 1)
x = pts[1][0]; y = pts[1][1]; z = pts[1][2]
x[0]= 6.; x[1] = 5.4; x[2]=4.8; x[3] = 2.5; x[4]  = 0.3
y[0]=10.; y[1]=0.036; y[2]=-5.;y[3]=0.21;y[4]=0.26;y[5]=7.
z[0]=1.; z[1]=1.; z[2]=1.;z[3]=1.;z[4]=1.;z[5]=1.; z[6]=1.

a = D.bezier( pts, 50 )
L = T.splitCurvatureRadius(a)
C.convertArrays2File([a]+L, 'out.plt')
