# - computeGrad2 (pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G
import KCore.test as test

# 3D STRUCT
ni = 30; nj = 40; nk = 2
t = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
C._initVars(t,'{centers:Density}=2*{centers:CoordinateX}+{centers:CoordinateX}*{centers:CoordinateY}')
P._computeGrad2(t, 'centers:Density')
test.testT(t, 1)

ni = 30; nj = 40; nk = 10
t = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
C._initVars(t,'{centers:Density}=2*{centers:CoordinateX}+{centers:CoordinateX}*{centers:CoordinateY}')
P._computeGrad2(t, 'centers:Density')
test.testT(t, 2)

# 3D NGON
ni = 30; nj = 40; nk = 2
t = G.cartNGon((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
C._initVars(t,'{centers:Density}=2*{centers:CoordinateX}+{centers:CoordinateX}*{centers:CoordinateY}')
P._computeGrad2(t, 'centers:Density')
test.testT(t, 3)

ni = 30; nj = 40; nk = 10
t = G.cartNGon((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
C._initVars(t,'{centers:Density}=2*{centers:CoordinateX}+{centers:CoordinateX}*{centers:CoordinateY}')
P._computeGrad2(t, 'centers:Density')
test.testT(t, 4)

# 2D STRUCT (nk=1)
ni = 30; nj = 40
t = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,1))
C._initVars(t,'{centers:Density}=2*{centers:CoordinateX}+{centers:CoordinateX}*{centers:CoordinateY}')
P._computeGrad2(t, 'centers:Density')
test.testT(t, 0)

# 2D NGON (z=0, nk=1)
ni = 30; nj = 40
t = G.cartNGon((0,0,0), (10./(ni-1),10./(nj-1),0), (ni,nj,1))
C._initVars(t,'{centers:Density}={centers:CoordinateX}')
P._computeGrad2(t, 'centers:Density')
C.convertPyTree2File(t, "out.cgns")
#test.testT(t, 5)
