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
t = C.newPyTree(['Base',1]); t[2][1][2].append(m)
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = P.computeGrad(t, 'Density')
t = C.initVars(t, 'centers:Pressure', F, ['gradxDensity','gradyDensity'])
t = P.computeNormGrad(t, 'centers:Pressure')
test.testT(t, 1)

#-----
# 2D
#-----
ni = 30; nj = 40
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,1))
m = C.initVars(m, 'Density', F, ['CoordinateX','CoordinateY'])
m = C.addBC2Zone(m,'ov','BCOverlap','imin')
t = C.newPyTree(['Base',2]); t[2][1][2].append(m)
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = P.computeGrad(t,'Density')
t = C.initVars(t, 'centers:Pressure', F, ['gradxDensity','gradyDensity'])
t = P.computeNormGrad(t,'centers:Pressure')
test.testT(t,2)

#-----
# 3D
#-----
ni = 30; nj = 40
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,2))
m = C.initVars(m, 'Density', F, ['CoordinateX','CoordinateY'])
t = C.newPyTree(['Base',3]); t[2][1][2].append(m)
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = C.fillEmptyBCWith(t, 'wall', 'BCWall', dim=2)
t = P.computeGrad(t, 'Density')
t = C.initVars(t, 'centers:Pressure', F, ['gradxDensity','gradyDensity'])
t = P.computeNormGrad(t, 'centers:Pressure')
test.testT(t,3)
