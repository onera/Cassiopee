# - streamSurf (pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G
import Geom.PyTree as D

ni = 30; nj = 40; nk = 5
m1 = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk)); m1[0] = 'cart1'
m2 = G.cart((5.5,0,0), (9./(ni-1),9./(nj-1),1), (ni,nj,nk)); m2[0] = 'cart2'
b = D.line((0.1,5.,0.1), (0.1,5.,3.9), N=10)
b = C.convertArray2Tetra(b)

t = C.newPyTree(['Base',m1,m2])
t = C.initVars(t, 'vx', 1.)
t = C.initVars(t, '{vy}=cos({CoordinateX})')
t = C.initVars(t, 'vz', 0.)
x0=0.1; y0=5.; z0=0.5
p = P.streamSurf(t, b, ['vx','vy','vz'])
C.convertPyTree2File(p, "out.cgns")
