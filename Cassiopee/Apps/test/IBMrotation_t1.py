# case - cylinder - IBM rotation - generate mesh and connectivity tree
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal
import Connector.PyTree as X
import Transform.PyTree as T
import Geom.PyTree as D
import Initiator.PyTree as I
import Generator.IBM as G_IBM
import Geom.IBM as D_IBM
import Connector.IBM as X_IBM
import RigidMotion.PyTree as R
import Connector.Mpi as Xmpi
import KCore.test as test
import math, os

LOCAL = test.getLocal()

## =======================
## =======================
TInf = 298.15
UInf = 0.0001
rad  = 0.0005
diam = rad*2
RoInf= 1.18
PInf = 101325

# necessaire pour les IBM
Model        = 'NSLaminar'
dimPb        = 2
frontType    = 1
vmin         = 11

# extension du maillage near-body jusqu au niveau de raffinement lvl_nearbody
dfar_nb  = 2*diam/2
dfar_ext = 5*diam/2
offset   = dfar_nb*0.25

rpm  = 100
OMG  = rpm*2*math.pi/60
LInf = diam*0.5
TInf = 298.15
aref = math.sqrt(1.4*287.058*298.15)
Vtip = OMG*LInf
Mtip = OMG*LInf/aref
Re   = RoInf*Vtip*diam/1.85e-05
## =======================
## =======================

##CREATING GEOMETRY
densityOfPts = 200
a = D.circle((0.,0.,0.), diam*0.5, N=densityOfPts)
  
C._initVars(a, '{CoordinateZ}=0')
tb           = C.newPyTree(["CYLINDER"]); tb[2][1][2] = Internal.getZones(a)
snear        = 5.0e-5
ModelTmp     = Model
D_IBM._setSnear(tb, snear)
D_IBM._setDfar(tb, dfar_ext)
D_IBM._setIBCType(tb,"noslip")
C._addState(tb, adim='dim1', UInf=UInf, TInf=TInf, PInf=PInf, LInf=diam,EquationDimension=dimPb, GoverningEquations=ModelTmp)
    
##CREATING OVERSET IBM MESH AROUND GEOMETRY
dfars=[]
for z in Internal.getZones(tb): dfars.append(dfar_nb)
t_ibm = G_IBM.generateIBMMesh(tb, vmin=vmin, dfars=dfars, dimPb=2, expand=2, ext=3)

R._setPrescribedMotion3(t_ibm ,'rot', axis_pnt=(0.,0.,0.), axis_vct=(0,0,1),omega=-OMG)
t_ibm[2][1][0]='CARTESIAN_NEARBODY'

t_out =None
tc_out=None
t_ibm, tc_ibm = X_IBM.prepareIBMData(tb, t_out, tc_out, t_ibm, frontType=frontType, cartesian=False)

C._rmBCOfType(t_ibm,'BCFarfield')
C._fillEmptyBCWith(t_ibm,'dummy','BCExtrapolate', dim=dimPb)

##GETTING EDGE OF OVERSET IBM MESH --> NEW GEOMETRY FOR BACKGROUND
ovs = C.extractBCOfType(t_ibm,'BCExtrapolate')

G._getVolumeMap(ovs)    
DZ = C.getMaxValue(ovs,'CoordinateZ')-C.getMinValue(ovs,'CoordinateZ')
C._initVars(ovs,'{centers:vol}={centers:vol}*%g'%(1./DZ))
ovs = T.subzone(ovs,(1,1,1),(-1,1,1))
C._initVars(ovs,'CoordinateZ',0)
ovs=C.convertArray2Hexa(ovs)
ovs=T.join(ovs)

for zs in Internal.getZones(ovs):
    snear = C.getMaxValue(zs,'centers:vol')
    snear = snear**(dimPb-1)
    D_IBM._setSnear(zs, snear)
    D_IBM._setDfar(zs, dfar_ext)
C._addState(ovs, adim='dim1', UInf=UInf, TInf=TInf, PInf=PInf, LInf=diam,EquationDimension=dimPb, GoverningEquations=Model)
tb_off = C.newPyTree(['Base', ovs])

##CREATING BACKGROUND CARTESIAN MESH
vmin    = 21
dfars=[]
for z in Internal.getZones(tb_off): dfars.append(dfar_ext)
t_off = G_IBM.generateIBMMesh(tb_off, vmin=vmin, dfars=dfars, dimPb=2, expand=2, ext=3)

t_off[2][1][0]='CARTESIAN_OFFBODY'
X._applyBCOverlaps(t_off,depth=2,loc='centers')
tc_off = C.node2Center(t_off)

cartesian=True
if Internal.getNodeFromType(t_off, "GridConnectivity1to1_t") is not None:
    Xmpi._setInterpData(t_off, tc_off, nature=1, loc='centers', storage='inverse', sameName=1, dim=dimPb, itype='abutting', order=2, cartesian=cartesian)
Xmpi._setInterpData(t_off, tc_off, nature=1, loc='centers', storage='inverse', sameName=1, sameBase=1, dim=dimPb, itype='chimera', order=2, cartesian=cartesian)
    

# ASSEMBLY
for a in [t_off, tc_off, t_ibm, tc_ibm]:
    C._addState(a, adim='dim1', UInf=UInf, TInf=TInf, PInf=PInf, LInf=diam,EquationDimension=dimPb, GoverningEquations=Model)
for t in [t_off,t_ibm]:
    C._initVars(t,"{centers:cellN#Static}={centers:cellN}")
    C._initVars(t,'centers:cellN',1.)
for var in ['cellN#Static',"cellN"]:
    C._cpVars(t_off,"centers:%s"%var,tc_off,var)
    C._cpVars(t_ibm,"centers:%s"%var,tc_ibm,var)

R._setPrescribedMotion3(t_ibm ,'rot', axis_pnt=(0.,0.,0.), axis_vct=(0,0,1),omega=-OMG)
R._setPrescribedMotion3(tc_ibm,'rot', axis_pnt=(0.,0.,0.), axis_vct=(0,0,1),omega=-OMG)


## BLANKING OF BACKGROUND MESH -- USING AN OFFSET OF IMMERSED GEOMETRY 
blankingBody= D.offsetSurface(Internal.getZones(tb_off), offset=-offset, pointsPerUnitLength=20000, algo=0, dim=dimPb)
tb_blank    = C.newPyTree(["BODY"]); tb_blank[2][1][2] = Internal.getZones(blankingBody)
R._setPrescribedMotion3(tb_blank ,'rot', axis_pnt=(0.,0.,0.), axis_vct=(0,0,1),omega=-OMG)
tb_blank = G.close(tb_blank,1.e-8)

if dimPb==2:
    zl = Internal.getZones(t_ibm)[0]
    DZ = C.getMaxValue(zl,'CoordinateZ')-C.getMinValue(zl,'CoordinateZ')
    T._addkplane(tb_blank)
    T._contract(tb_blank,  (0,0,0), (1,0,0), (0,1,0), DZ)
R._copyGrid2GridInit(tb_blank)    
C.convertPyTree2File(tb_blank, LOCAL+'/bodiesBlank.cgns')

# suppress static BCOverlap in t_off
C._rmBCOfType(t_off, 'BCFarfield')
C._fillEmptyBCWith(t_off, 'nref', 'BCFarfield', dim=dimPb)
C._rmBCOfType(t_off, 'BCOverlap')
# t_ibm
C._fillEmptyBCWith(t_ibm,'dummy','BCExtrapolate',dim=dimPb)# corners
C._rmBCOfType(t_ibm, 'BCOverlap')
C._fillEmptyBCWith(t_ibm, 'dummy', 'DUMMY', dim=dimPb) # overlap static
C._rmBCOfType(t_ibm, 'BCExtrapolate')
C._fillEmptyBCWith(t_ibm, 'ovst_dyn', 'BCOverlap', dim=dimPb)
C._rmBCOfType(t_ibm, 'DUMMY')

# Init field
vars = ['Density', 'MomentumX', 'MomentumY', 'MomentumZ', 'EnergyStagnationDensity']
C._initVars(t_off, "centers:TurbulentDistance", 1e3)

t  = C.mergeTrees(t_ibm, t_off)
tc = C.mergeTrees(tc_ibm, tc_off)

for b in Internal.getBases(t):
    state = Internal.getNodeFromType(b, 'ReferenceState_t')
    #Model = Internal.getNodeFromName(b, 'GoverningEquations')
    #Model = Internal.getValue(Model)
    if Model == 'NSTurbulent': allvars = vars + ['TurbulentSANuTildeDensity']
    else: allvars = vars
    for v in allvars:
        node = Internal.getNodeFromName(state, v)
        if node is not None:
            val = float(node[1][0])
            C._initVars(b, 'centers:'+v, val)
            
R._copyGrid2GridInit(t)
C.convertPyTree2File(t, LOCAL+'/t.cgns')
R._copyGrid2GridInit(tc)
C.convertPyTree2File(tc, LOCAL+'/tc.cgns')

Internal._rmNodesByName(t, '.Solver#Param')
Internal._rmNodesByName(t, '.Solver#ownData')
Internal._rmNodesByName(tc, '.Solver#Param')
Internal._rmNodesByName(tc, '.Solver#ownData')

test.testT(t, 1)
test.testT(tc, 2)

