# - computeVariables (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
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
t = C.newPyTree(['Base',1]); t[2][1][2].append(m)
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = P.computeVariables(t, ['Mach', 'Pressure'])
test.testT(t,1)
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
m = C.addBC2Zone(m,'ov','BCOverlap','imin')
m = C.addBC2Zone(m,'wall','BCWall','imax')
t = C.newPyTree(['Base',2]); t[2][1][2].append(m)
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = P.computeVariables(t, ['Mach', 'Pressure'])
test.testT(t,2)

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
m = C.addBC2Zone(m,'ov','BCOverlap','imin')
m = C.addBC2Zone(m,'wall','BCWall','imax')
t = C.newPyTree(['Base']); t[2][1][2].append(m)
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = P.computeVariables(t, ['Mach', 'Pressure'])
test.testT(t,3)
#
# calcul des variables aux centres
#
ni = 30; nj = 40
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,2))
m = C.initVars(m,'Density',initDensity,['CoordinateX','CoordinateY','CoordinateZ'])
m = C.initVars(m,'MomentumX',initMomentum,['CoordinateX','CoordinateY','CoordinateZ'])
m = C.initVars(m,'MomentumY',initMomentum,['CoordinateX','CoordinateY','CoordinateZ'])
m = C.initVars(m,'MomentumZ',initMomentum,['CoordinateX','CoordinateY','CoordinateZ'] )
m = C.initVars(m,'EnergyStagnationDensity',initEnergy,['CoordinateX','CoordinateY','CoordinateZ'])
m = C.node2Center(m, Internal.__FlowSolutionNodes__)
m = C.addBC2Zone(m,'ov','BCOverlap','imin')
m = C.addBC2Zone(m,'wall','BCWall','imax')
t = C.newPyTree(['Base']); t[2][1][2].append(m)
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = P.computeVariables(t, ['centers:Mach', 'centers:Pressure'])
test.testT(t,4)
