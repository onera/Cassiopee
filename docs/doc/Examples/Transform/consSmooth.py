# - consSmooth (array) -
import Transform as T
import Converter as C
import Generator as G
import Geom as D

# 2D : square
l1 = D.line((0.,0.,0.), (0.,1.,0.), N=5)
l2 = D.line((0.,1.,0.), (1.,1.,0.), N=5)
l3 = D.line((1.,1.,0.), (1.,0.,0.), N=5)
l4 = D.line((1.,0.,0.), (0.,0.,0.), N=5)

a = T.join([l1,l2,l3,l4])
b = T.consSmooth(a, sweeps=5)
C.convertArrays2File([a,b], "out1.plt")

# 3D : cube
e = D.box((0,0,0), (10,10,10), N=10, ntype='QUAD')
e = C.convertArray2Tetra(e, split='withBarycenters')
e = G.close(e)
e = T.consSmooth(e, sweeps=5, omega=0.1)
C.convertArrays2File(e, "out2.plt")
