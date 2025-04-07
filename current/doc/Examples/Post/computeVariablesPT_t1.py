# - computeVariables (pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G
import KCore.test as test

def initDensity(x,y,z): return 1.
def initMomentum(x,y,z): return 0.1+3*x + 5*y +2*z
def initEnergy(x,y,z): return 1e5+1.5*x*x + y*y + 3.*z*z

#-----
# 1D
#-----
ni = 30
m = G.cart((0,0,0), (10./(ni-1),1,1), (ni,1,1))
m = C.initVars(m,'Density',initDensity,['CoordinateX','CoordinateY','CoordinateZ'])
m = C.initVars(m,'centers:Density',1.)
m = C.initVars(m,'MomentumX',initMomentum,['CoordinateX','CoordinateY','CoordinateZ'])
m = C.initVars(m,'MomentumY',0.)
m = C.initVars(m,'MomentumZ',0.)
m = C.initVars(m,'EnergyStagnationDensity',initEnergy,['CoordinateX','CoordinateY','CoordinateZ'])
m = P.computeVariables(m, ['Pressure','Mach'])
test.testT(m,1)

#-----
# 2D
#-----
ni = 30; nj = 40
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,1))
m = C.initVars(m,'Density',initDensity,['CoordinateX','CoordinateY','CoordinateZ'])
m = C.initVars(m,'centers:Density',1.)
m = C.initVars(m,'MomentumX',initMomentum,['CoordinateX','CoordinateY','CoordinateZ'])
m = C.initVars(m,'MomentumY',initMomentum,['CoordinateX','CoordinateY','CoordinateZ'])
m = C.initVars(m,'MomentumZ',0.)
m = C.initVars(m,'EnergyStagnationDensity',initEnergy,['CoordinateX','CoordinateY','CoordinateZ'])
m = P.computeVariables(m, ['Mach', 'Pressure'])
test.testT(m,2)

#-----
# 3D
#-----
ni = 30; nj = 40
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,2))
m = C.initVars(m,'centers:Density',1.)
m = C.initVars(m,'Density',initDensity,['CoordinateX','CoordinateY','CoordinateZ'])
m = C.initVars(m,'MomentumX',initMomentum,['CoordinateX','CoordinateY','CoordinateZ'])
m = C.initVars(m,'MomentumY',initMomentum,['CoordinateX','CoordinateY','CoordinateZ'])
m = C.initVars(m,'MomentumZ',initMomentum,['CoordinateX','CoordinateY','CoordinateZ'] )
m = C.initVars(m,'EnergyStagnationDensity',initEnergy,['CoordinateX','CoordinateY','CoordinateZ'])
m = P.computeVariables(m, ['Mach', 'Pressure'])
test.testT(m,3)
