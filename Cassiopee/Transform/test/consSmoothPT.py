# - consSmooth (pyTree) -
import Transform.PyTree as T
import Converter.PyTree as C
import Geom.PyTree as D
import Generator.PyTree as G
import Converter.Internal as I

# 2D : square
l1 = D.line((0.,0.,0.), (0.,1.,0.), N=5)
l2 = D.line((0.,1.,0.), (1.,1.,0.), N=5)
l3 = D.line((1.,1.,0.), (1.,0.,0.), N=5)
l4 = D.line((1.,0.,0.), (0.,0.,0.), N=5)

a = T.join([l1,l2,l3,l4]); I.setName(a, 'carreOriginal')
b = T.consSmooth(a, sweeps=5); I.setName(b, 'carreLisse')
C.convertPyTree2File([a,b], "out1.cgns")

# 3D : cube
e = D.box((0.0, 0.0, 0.0), (10, 10, 10), N=5, ntype='QUAD')
e = C.convertArray2Tetra(e, split='withBarycenters')
e = G.close(e)
e = T.consSmooth(e, sweeps=5, omega=0.1)
C.convertPyTree2File(e, "out2.cgns")
