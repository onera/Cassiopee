# - nurbs (array) -
import Geom as D
import Converter as C
import Generator as G

a = D.polyline ([(4.1,0.1,1.1), (1.1,0.2,1.2), (1.1,1.3,1.3),
                 (1.1,1.5,1.4), (4.5,2.5,1.5), (5.6,1.5,1.6),
                 (6.7,1.7,1.7), (7.8,0.8,1.8), (8.9,-1.9,1.9), (9,0,1)])
a = C.initVars(a,'W',1.)
C.convertArrays2File([a],'in.plt')
b = D.nurbs(a,"W", 4, N=100)
c = D.nurbs(a,"W", 4, density=10.)
C.convertArrays2File([b,c], 'out.plt')

ni = 10; nj = 10
a = G.cart((0,0,0), (1,1,1), (ni,nj,1))
C.setValue(a, (1,1,1), [1.,1.,1.])
C.setValue(a, (1,2,1), [1.,2.,1.])
C.setValue(a, (1,3,1), [1.,3.,1.])
C.setValue(a, (1,4,1), [1.,4.,1.])
C.setValue(a, (2,1,1), [2.,1.,2.])
C.setValue(a, (2,2,1), [2.,2.,5.])
C.setValue(a, (2,3,1), [2.,3.,5.])
C.setValue(a, (2,4,1), [2.,4.,2.])
C.setValue(a, (3,1,1), [3.,1.,2.])
C.setValue(a, (3,2,1), [3.,2.,5.])
C.setValue(a, (3,3,1), [3.,3.,12.])
C.setValue(a, (3,4,1), [3.,4.,2.])
C.setValue(a, (4,1,1), [4.,1.,2.])
C.setValue(a, (4,2,1), [4.,2.,5.])
C.setValue(a, (4,3,1), [4.,3.,5.])
C.setValue(a, (4,4,1), [4.,4.,2.])
C.setValue(a, (6,8,1), [4.,6.,14.])
C.setValue(a, (8,6,1), [4.,6.,-4.])
a = C.initVars(a,"W",1.)
a[1][3,6]=7; a[1][3,14]=9.
d = D.nurbs(a, "W", 4, N=100, M=100)
e = D.nurbs(a, "W", 4, density=20.)
C.convertArrays2File([a],'in2.plt')
C.convertArrays2File([d,e],'out2.plt')
