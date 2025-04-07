# - extractBCFields (pyTree) -
import Converter.Internal as Internal
import Converter.PyTree as C
import Generator.PyTree as G

ni = 10; nj = 10; nk = 10
nfaces = (nj-1)*(nk-1)
a = G.cart((0,0,0), (1,1,1), (ni,nj,nk))
C._addBC2Zone(a, 'wall', 'BCWall', 'imin')

C._initVars(a,'centers:VelocityX={centers:CoordinateX}')
C._initVars(a,'centers:Density=1.05')

b = Internal.getNodeFromName2(a, 'wall')
d = Internal.newBCDataSet(name='BCDataSet', value='UserDefined',
                          gridLocation='FaceCenter', parent=b)
d1 = Internal.newBCData('BCNeumann', parent=d)
d = Internal.newDataArray('Density', value=nfaces*[1.], parent=d1)
d = Internal.newDataArray('MomentumX', value=nfaces*[0.3], parent=d1)

# Get data array node list
#varList=['Density','MomentumX']
varList=None
res = C.extractBCFields(a,varList=varList)
print('variables = ', res[0])
print('fields = ',  res[1])
print('indices = ', res[2])
