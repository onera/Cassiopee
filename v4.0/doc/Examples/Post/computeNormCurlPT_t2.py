# - computeNormCurl(pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G
import KCore.test as test

def F(x,y,z): return 12*y*y + 4

#-----
# 2D
#-----
ni = 30; nj = 40
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,1))
m = C.initVars(m,'F1',F,['CoordinateX','CoordinateY','CoordinateZ'])
m = C.addVars(m,'F2'); m = C.addVars(m,'F3')
varname = ['F1','F2','F3']
m = P.computeCurl(m, varname)
m = C.initVars(m,'centers:F4',F,['centers:rotx','centers:roty','centers:rotz'])
m = C.addVars(m,'centers:F5'); m = C.addVars(m,'centers:F6')
m = C.addBC2Zone(m,'ov','BCOverlap','imin')
t = C.newPyTree(['Base',2]); t[2][1][2].append(m)
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
varname = ['centers:F4','centers:F5','centers:F6']
t = P.computeNormCurl(t, varname)
test.testT(t,2)

#-----
# 3D
#-----
ni = 30; nj = 40
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,2))
m = C.initVars(m,'F1',F,['CoordinateX','CoordinateY','CoordinateZ'])
m = C.addVars(m,'F2'); m = C.addVars(m,'F3')
varname = ['F1','F2','F3']
m = P.computeCurl(m, varname)
m = C.initVars(m,'centers:F4',F,['centers:rotx','centers:roty','centers:rotz'])
m = C.addVars(m,'centers:F5'); m = C.addVars(m,'centers:F6')
m = C.addBC2Zone(m,'ov','BCOverlap','imin')
t = C.newPyTree(['Base',3]); t[2][1][2].append(m)
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = C.fillEmptyBCWith(t,'wall','BCWall', dim=2)
varname = ['centers:F4','centers:F5','centers:F6']
t = P.computeNormCurl(t, varname)
test.testT(t,3)
