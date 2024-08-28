# - map (array) -
import Geom as D
import Generator as G
import Converter as C

# Map on a curve
l = D.line( (0,0,0), (1,1,0) )
Ni = 10
d = G.cart( (0,0,0), (1./(Ni-1),1.,1.), (Ni,1,1) )
m = G.map(l, d)
C.convertArrays2File([m], "out1.plt")

# Map on a structured surface
ni = 2; nj = 3
a = G.cart((0,0,0), (1,1,1), (ni,nj,1))
C.setValue(a, (1,1,1), [1.,1.,2.])
C.setValue(a, (1,2,1), [1.,2.,5.])
C.setValue(a, (1,3,1), [1.,3.,2.])
C.setValue(a, (2,1,1), [2.,1.,2.])
C.setValue(a, (2,2,1), [2.,2.,5.])
C.setValue(a, (2,3,1), [2.,3.,2.])
b = D.bezier(a, 10, 10)
Ni = 50; Nj = 30
d = G.cart( (0,0,0), (1./(Ni-1),1./(Nj-1),1.), (Ni,Nj,1) )
d = G.enforceX(d, 0.5, 0.01, (10,20))
d = G.enforceY(d, 0.5, 0.01, (10,20))
b = G.map(b, d)
C.convertArrays2File(b, "out2.plt")

# Map in a direction
a = G.cylinder((0,0,0), 0.5, 2., 0, 60, 1., (20,20,1))
Ni = 10
d = G.cart( (0,0,0), (1./(Ni-1),1.,1.), (Ni,1,1) )
d = G.enforcePlusX(d, 0.01, (10,20))
a = G.map(a, d, 2)
C.convertArrays2File(a, "out3.plt")
