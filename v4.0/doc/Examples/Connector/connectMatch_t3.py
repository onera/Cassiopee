# connectMatch sur meme bloc (array)
import Generator as G
import Connector as X
import Converter as C
import Geom as D
import Transform as T
import KCore.test as test

# 3D raccord i=1 avec i = imax
a = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0., 10., (50,50,30))
res = X.connectMatch(a,a,sameZone=1)
test.testO(res,1)

# 3D raccord i = 1 partiel profil NACA
msh = D.naca(12., 5001)
msh2 = D.line((1.,0.,0.),(2.,0.,0.),5001); msh = T.join(msh, msh2)
msh2 = D.line((2.,0.,0.),(1.,0.,0.),5001); msh = T.join(msh2, msh)
Ni = 300; Nj = 50
distrib = G.cart((0,0,0), (1./(Ni-1), 0.5/(Nj-1),1), (Ni,Nj,1))
naca = G.hyper2D(msh, distrib, "C")
a = T.addkplane(naca)
res = X.connectMatch(a,a,sameZone=1)
test.testO(res,2)

# 2D raccord i=1 avec i = imax
a = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0., 10., (50,50,1))
ca = C.array('cellN',49,49,0)
ca = C.initVars(ca, 'cellN', 1.)
a = C.addVars(a, 'cellN')
res = X.connectMatch(a,a,sameZone=1,dim=2)
test.testO(res,3)

# 2D raccord i = 1 partiel
res = X.connectMatch(naca,naca,sameZone=1,dim=2)
test.testO(res,4)
