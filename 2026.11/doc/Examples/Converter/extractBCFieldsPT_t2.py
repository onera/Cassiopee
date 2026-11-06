# - extractBCFields (pyTree) -
import Converter.Internal as Internal
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

ni = 10; nj = 10; nk = 10
nfaces = (nj-1)*(nk-1)
a = G.cart((0,0,0), (1,1,1), (ni,nj,nk))
C._addBC2Zone(a, 'wall', 'BCWall', 'imin')

C._initVars(a, 'centers:VelocityX={centers:CoordinateX}')
C._initVars(a, 'centers:Density=1.05')
a=C.convertArray2NGon(a)
b = Internal.getNodeFromName2(a, 'wall.1')
d = Internal.newBCDataSet(name='BCDataSet', value='UserDefined',
                          gridLocation='FaceCenter', parent=b)
d1 = Internal.newBCData('BCNeumann', parent=d)
d = Internal.newDataArray('Density', value=nfaces*[1.], parent=d1)
d = Internal.newDataArray('MomentumX', value=nfaces*[0.3], parent=d1)

varList=None
res = C.extractBCFields(a, varList=varList)
test.testO(res[0],11)
test.testO(res[1],12)
test.testO(res[2],13)

varList=['Density','MomentumX']
res = C.extractBCFields(a,varList=varList)
test.testO(res[0],21)
test.testO(res[1],22)
test.testO(res[2],23)

a = C.convertArray2NGon(a)
varList=None
res = C.extractBCFields(a,varList=varList)
test.testO(res[0],31)
test.testO(res[1],32)
test.testO(res[2],33)

varList=['Density','MomentumX']
res = C.extractBCFields(a,varList=varList)
test.testO(res[0],41)
test.testO(res[1],42)
test.testO(res[2],43)
