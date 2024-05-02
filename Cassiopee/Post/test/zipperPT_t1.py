# - zipper (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T
import KCore.test as test

# cylindre
ni = 30; nj = 40; nk = 1
m1 = G.cylinder((0,0,0), 1, 5, 0., 360., 10., (50,50,1))
m1 = C.initVars(m1, 'cellN', 1.)
m1 = C.initVars(m1, 'centers:ichim', 1.)

# Set cellN = 2 (interpolated points) to boundary
s = T.subzone(m1, (1,nj,1),(ni,nj,nk))
s = C.initVars(s, 'cellN', 2)
m1 = T.patch(m1, s, (1,nj,1))

ni1 =  m1[1][0][0]; nj1 =  m1[1][0][1]; nk1 = nk
s = T.subzone(m1, (1,1,1),(ni1,1,nk1))
s = C.initVars(s, 'cellN', 2)
m1 = T.patch(m1, s, (1,1,1))

# carre
m2 = G.cart((0,0,0), (10./(ni-1),10./(nj-1),-1), (ni,nj,1))
m2 = C.initVars(m2, 'cellN',1.); m2 = C.initVars(m2, 'Density', 1.2)
test.testT(m2)
