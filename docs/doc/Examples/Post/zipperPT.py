# - zipper (pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G
import Transform.PyTree as T

# cylindre
ni = 30; nj = 40; nk = 1
m1 = G.cylinder((0,0,0), 1, 5, 0., 360., 10., (ni,nj,nk))
m1 = C.addVars(m1, 'Density'); m1 = C.initVars(m1,'cellN',1)

# Set cellN = 2 (interpolated points) to boundary
s = T.subzone(m1, (1,nj,1),(ni,nj,nk))
s = C.initVars(s, 'cellN', 2)
m1 = T.patch(m1, s, (1,nj,1))
s = T.subzone(m1, (1,1,1),(ni,1,nk))
s = C.initVars(s, 'cellN', 2)
m1 = T.patch(m1, s, (1,1,1))

# carre
ni = 30; nj = 40
m2 = G.cart((0,0,0), (10./(ni-1),10./(nj-1),-1), (ni,nj,1))
m2 = C.initVars(m2, 'Density', 1.2); m2 = C.initVars(m2, 'cellN', 1.)

t = C.newPyTree(['Base',2]); t[2][1][2] += [m1, m2]
z = P.zipper(t); z[0] = 'zipper'; t[2][1][2].append(z)
C.convertPyTree2File(t, 'out.cgns')
