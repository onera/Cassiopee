# - computeDiv (pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G
import KCore.test as test

def Fx(x,y,z): return 2*x+x*y
def Fy(x,y,z): return 4.*y
def Fz(x,y,z): return x*y+z*z

#-----
# 2D
#-----
ni = 30; nj = 40
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,1))
m = C.initVars(m, 'fldX', Fx, ['CoordinateX','CoordinateY','CoordinateZ'])
m = C.initVars(m, 'fldY', Fy, ['CoordinateX','CoordinateY','CoordinateZ'])
m = C.initVars(m, 'fldZ', Fz, ['CoordinateX','CoordinateY','CoordinateZ'])
m = C.addBC2Zone(m,'ov','BCOverlap','imin')
t = C.newPyTree(['Base',2]); t[2][1][2].append(m)
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = P.computeDiv(t, 'fld')
t = C.initVars(t, 'centers:tmpX', Fx, ['centers:CoordinateX','centers:CoordinateY','centers:divfld'])
t = C.initVars(t, 'centers:tmpY', Fy, ['centers:CoordinateX','centers:CoordinateY','centers:divfld'])
t = C.initVars(t, 'centers:tmpZ', Fz, ['centers:CoordinateX','centers:CoordinateY','centers:divfld'])
t = P.computeDiv(t,'centers:tmp')
test.testT(t,1)

#-----
# 3D
#-----
ni = 30; nj = 40
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,2))
m = C.initVars(m, 'fldX', Fx, ['CoordinateX','CoordinateY','CoordinateZ'])
m = C.initVars(m, 'fldY', Fy, ['CoordinateX','CoordinateY','CoordinateZ'])
m = C.initVars(m, 'fldZ', Fz, ['CoordinateX','CoordinateY','CoordinateZ'])
t = C.newPyTree(['Base',3]); t[2][1][2].append(m)
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = C.fillEmptyBCWith(t, 'wall', 'BCWall', dim=2)
t = P.computeDiv(t, 'fld')
t = C.initVars(t, 'centers:tmpX', Fx, ['centers:CoordinateX','centers:CoordinateY','centers:divfld'])
t = C.initVars(t, 'centers:tmpY', Fy, ['centers:CoordinateX','centers:CoordinateY','centers:divfld'])
t = C.initVars(t, 'centers:tmpZ', Fz, ['centers:CoordinateX','centers:CoordinateY','centers:divfld'])
t = P.computeDiv(t,'centers:tmp')
test.testT(t,2)
