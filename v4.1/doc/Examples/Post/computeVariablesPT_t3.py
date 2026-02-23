# - computeVariables (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Post.PyTree as P
import Generator.PyTree as G
import KCore.test as test

#-------------------------
# Test reference state
#-------------------------
ni = 31; dh = 10./(ni-1)
m = G.cart((0,0,0), (dh,dh,1), (ni,ni,1))
m = C.initVars(m,'Density',1.)
m = C.initVars(m,'{MomentumX}=0.1+3.*{CoordinateX}+5.*{CoordinateY}+2.*{CoordinateZ}')
m = C.initVars(m,'MomentumY',0.)
m = C.initVars(m,'MomentumZ',0.)
m = C.initVars(m,'{EnergyStagnationDensity}=1.e5+1.5*{CoordinateX}**2+{CoordinateY}**2+3.*{CoordinateZ}**2')
t = C.newPyTree(['Base',1]); t[2][1][2].append(m)
t = C.addState(t, adim='adim1', MInf=0.6)
t = P.computeVariables(t, ['Mach', 'Pressure'])
test.testT(t,1)
