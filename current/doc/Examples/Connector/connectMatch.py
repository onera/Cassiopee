# - connectMatch (array) -
import Generator as G
import Connector as X
import Geom as D
import Transform as T
import Converter as C
# 3D raccord i = 1 partiel profil NACA
msh = D.naca(12., 5001)
msh2 = D.line((1.,0.,0.),(2.,0.,0.),5001); msh = T.join(msh, msh2)
msh2 = D.line((2.,0.,0.),(1.,0.,0.),5001); msh = T.join(msh2, msh)
Ni = 300; Nj = 50
distrib = G.cart((0,0,0), (1./(Ni-1), 0.5/(Nj-1),1), (Ni,Nj,1))
naca = G.hyper2D(msh, distrib, "C")
res = X.connectMatch(naca,naca,sameZone=1,dim=2)
C.convertArrays2File([naca],"out.plt")
print(res)
