# case - cylinder - automatic rectilinear mesh generator
import Connector.IBM as X_IBM
import Converter.Internal as Internal
import Converter.PyTree as C
import Generator.PyTree as G
import Geom.IBM as D_IBM
import Geom.PyTree as D
import KCore.test as test
import Post.PyTree as P
import os

LOCAL = test.getLocal()

## =======================
## =======================
TInf = 298.15
UInf = 0.0001
rad  = 0.5
diam = rad*2
RoInf= 1.18
PInf = 101325
LInf = diam*0.5
TInf = 298.15

Model        = 'NSLaminar'
dimPb        = 2
frontType    = 1
vmin         = 11
dfar         = 2*diam/2
## =======================
## =======================

##CREATING TB GEOMETRY
densityOfPts = 200
a = D.circle((0.,0.,0.), diam*0.5, N=densityOfPts)

C._initVars(a, '{CoordinateZ}=0')
tb           = C.newPyTree(["CYLINDER"]); tb[2][1][2] = Internal.getZones(a)
snear        = 1.0e-2
ModelTmp     = Model
D_IBM._setSnear(tb, snear)
D_IBM._setDfar(tb, dfar)
D_IBM._setIBCType(tb,"noslip")
C._addState(tb, adim='dim1', UInf=UInf, TInf=TInf, PInf=PInf, LInf=diam,EquationDimension=dimPb, GoverningEquations=ModelTmp)

##CREATING TB RECTILINEAR REGION
Lx    = dfar*3
hlocal= Lx/10
tbOneOver  = G.cart((-dfar*1.5,-dfar*1.5,0.), (hlocal,hlocal,0.1), (10,10,1))
tbOneOver2 = P.exteriorFaces(tbOneOver)
tbOneOver  = C.newPyTree(["Base_OneOver"]); tbOneOver[2][1][2] = Internal.getZones(tbOneOver2)
C.convertPyTree2File(tbOneOver, LOCAL+'/tbOneOver.cgns')

##ADD DIRECTION OF ONE OVER IN TB RECTILINEAR REGION
FileNameOneOver = LOCAL+'/tbOneOver.cgns'
listOneOver     = [[2,1,1]]
X_IBM._addOneOverLocally(FileNameOneOver, listOneOver)
tbOneOver = C.convertFile2PyTree(FileNameOneOver)

##IBM PREP
t,tc=X_IBM.prepareIBMData(tb         , None      , None     , tbox=tbOneOver,
                          snears=0.01, dfars=dfar, vmin=vmin, cartesian=False)
os.remove(FileNameOneOver)

####
# The following lines are to avoid regression since the bug fix for duplicate information in tc
####
for b in Internal.getBases(tc):
    for z in Internal.getZones(b):
        pos = 0
        z2 = Internal.copyRef(z)
        for zs in z2[2]:
            if 'ID' in zs[0] or 'IBCD' in zs[0]:
                Internal.addChild(z, zs, pos)
                pos +=2
            else:
                pos += 1
####

##NON-REGRESSION CHECK
Internal._rmNodesByName(t, '.Solver#Param')
Internal._rmNodesByName(t, '.Solver#ownData')
Internal._rmNodesByName(tc, '.Solver#Param')
Internal._rmNodesByName(tc, '.Solver#ownData')

test.testT(t, 1)
test.testT(tc, 2)

#C.convertPyTree2File(t ,LOCAL+'/t_check.cgns')
#C.convertPyTree2File(tc,LOCAL+'/tc_check.cgns')
