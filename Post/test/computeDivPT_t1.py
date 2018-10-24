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
m = P.computeDiv(m, 'fld')
m = C.initVars(m, 'centers:tmpX', Fx, ['centers:CoordinateX','centers:CoordinateY','centers:divfld'])
m = C.initVars(m, 'centers:tmpY', Fy, ['centers:CoordinateX','centers:CoordinateY','centers:divfld'])
m = C.initVars(m, 'centers:tmpZ', Fz, ['centers:CoordinateX','centers:CoordinateY','centers:divfld'])
m = P.computeDiv(m, 'centers:tmp')
test.testT(m,1)

#-----
# 3D
#-----
ni = 30; nj = 40
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,2))
m = C.initVars(m, 'fldX', Fx, ['CoordinateX','CoordinateY','CoordinateZ'])
m = C.initVars(m, 'fldY', Fy, ['CoordinateX','CoordinateY','CoordinateZ'])
m = C.initVars(m, 'fldZ', Fz, ['CoordinateX','CoordinateY','CoordinateZ'])
m = P.computeDiv(m, 'fld')
m = C.initVars(m, 'centers:tmpX', Fx, ['centers:CoordinateX','centers:CoordinateY','centers:divfld'])
m = C.initVars(m, 'centers:tmpY', Fy, ['centers:CoordinateX','centers:CoordinateY','centers:divfld'])
m = C.initVars(m, 'centers:tmpZ', Fz, ['centers:CoordinateX','centers:CoordinateY','centers:divfld'])
m = P.computeDiv(m, 'centers:tmp')
test.testT(m,2)
