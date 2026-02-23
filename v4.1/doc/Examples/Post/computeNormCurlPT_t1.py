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
varname = ['centers:F4','centers:F5','centers:F6']
m = P.computeNormCurl(m, varname)
test.testT(m,2)

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
varname = ['centers:F4','centers:F5','centers:F6']
m = P.computeNormCurl(m, varname)
test.testT(m,3)
