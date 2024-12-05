# - extractionIBM a la paroi (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Connector.PyTree as X
import Geom.PyTree as D
import Post.PyTree as P
import Dist2Walls.PyTree as DTW
import Transform.PyTree as T
import Initiator.PyTree as I
import Converter.Internal as Internal
import Connector.ToolboxIBM as IBM
import KCore.test as test

N = 41
a = G.cart((0,0,0),(1./(N-1),1./(N-1),1./(N-1)),(N,N,N))
xm = 0.5*N/(N-1)
s = D.sphere((xm,xm,xm),0.1,N=20)
s = C.convertArray2Tetra(s); s = G.close(s)
s = T.splitNParts(s,4)
tb = C.newPyTree(['SPHERE']); tb[2][1][2]=s
C._addState(tb, 'EquationDimension',3)
C._addState(tb, 'GoverningEquations', 'NSTurbulent')

famNames=[]
for base in Internal.getBases(tb):
    baseName = Internal.getName(base)
    i=1
    for s in Internal.getZones(base):
        famName=baseName+'_'+str(i)
        famNames.append(famName)
        Internal._createChild(base,famName, 'Family_t',None)
        C._tagWithFamily(s,famName)
        i+=1

t = C.newPyTree(['Base', a])

# Dist2Walls
DTW._distance2Walls(t, tb, type='ortho', loc='centers')
# Gradient de distance localise en centres => normales
t = P.computeGrad(t, 'centers:TurbulentDistance')
I._initConst(t, MInf=0.2, loc='centers')
t,tc=IBM.prepareIBMData_legacy(t, tb, DEPTH=2, frontType=1)
z = IBM.extractIBMWallFields(tc, tb=tb, famZones=[famNames[0]])
test.testT(z,1)
