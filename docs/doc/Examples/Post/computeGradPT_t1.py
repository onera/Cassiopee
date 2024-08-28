# - computeGrad (pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G
import KCore.test as test

def F(x,y): return 2*x+x*y

#-----
# 1D
#-----
ni = 30
m = G.cart((0,0,0), (10./(ni-1),1,1), (ni,1,1))
m = C.initVars(m, 'Density', F, ['CoordinateX','CoordinateY'])
m = P.computeGrad(m,'Density')
m = C.initVars(m, 'centers:Pressure', F, ['centers:gradxDensity','centers:gradyDensity'])
m = P.computeGrad(m, 'centers:Pressure')
test.testT(m,1)

#-----
# 2D
#-----
ni = 30; nj = 40
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,1))
m = C.initVars(m, 'Density', F, ['CoordinateX', 'CoordinateY'])
m = P.computeGrad(m, 'Density')
m = C.initVars(m, 'centers:Pressure', F, ['centers:gradxDensity','centers:gradyDensity'])
m = P.computeGrad(m, 'centers:Pressure')
test.testT(m,2)

#-----
# 3D
#-----
ni = 30; nj = 40
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,2))
m = C.initVars(m, 'Density', F, ['CoordinateX','CoordinateY'])
m = P.computeGrad(m, 'Density')
m = C.initVars(m, 'centers:Pressure', F, ['centers:gradxDensity','centers:gradyDensity'])
m = P.computeGrad(m, 'centers:Pressure')
test.testT(m,3)
