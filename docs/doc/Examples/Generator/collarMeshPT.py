# - collarMesh (pyTree)
import Converter.PyTree as C
import Geom.PyTree as D
import Transform.PyTree as T
import Generator.PyTree as G

s1 = D.sphere((0.,0.,0.),1,20)
s2 = T.translate(s1,(1.2,0.,0.)); s2 = T.homothety(s2,(0,0,0),0.5)
dhj = G.cart((0.,0.,0.),(1.e-2,1,1),(21,1,1))
dhk = G.cart((0.,0.,0.),(1.e-2,1,1),(11,1,1))
a = G.collarMesh(s1,s2, dhj,dhk,niterj=100,niterk=100,type='union')
C.convertPyTree2File(a,"out.cgns")
