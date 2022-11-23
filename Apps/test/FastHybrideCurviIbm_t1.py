# - Fast.IBM -
# NS, para, frontType=2
import Apps.Fast.IBM as App
import Fast.PyTree as Fast
import FastS.PyTree as FastS
import Converter.PyTree as C
import Initiator.PyTree as Initiator
import Generator.PyTree as G
import Connector.PyTree as X
import Converter.Internal as Internal
import math
import numpy as np
import KCore.test as test
test.TOLERANCE = 1.e-6

LOCAL = test.getLocal()

#load case 2d pour maillage octree 
tb_2d = C.convertFile2PyTree('cavityCase.cgns')

# Prepare cas 2d IBM
t_2d, tc_2d = App.prepare1_IM('cavityCase.cgns', t_out=LOCAL+'/t.cgns', tc_out=LOCAL+'/tc.cgns', frontType=42 ,cleanCellN=False)

#extrusion 3D IBM
extrusion ='cart'; span = 0.01
t_3d, tb_3d = App.extrudeCartesian(t_2d, tb_2d, extrusion=extrusion, NPas=5, span=span)

#calcul interp IBM 3D
interpDataType = 1 # on suppose maillage non cartesion pour interp
order          = 2
t_3d, tc_3d = App.prepare1_IM(tb_3d,None, None, t_in=t_3d, extrusion=extrusion, interpDataType=interpDataType, order=order)

Internal._rmNodesFromType(tc_2d,'Rind_t')
Internal._rmNodesFromName(tc_2d,Internal.__GridCoordinates__)
test.testT(tc_2d, 1)

#Maillage plaque plane Curviligne sans ghost
a = G.cart((0,0,0), (0.005,0.005,0.01), (200,100,5))
C._addBC2Zone(a, 'farfield', 'BCFarfield', 'jmax')
C._addBC2Zone(a, 'inflow', 'BCInflow',     'imin')
C._addBC2Zone(a, 'outflow','BCOutflow',    'imax')
C._addBC2Zone(a, 'wallIn','BCWallViscous',  [1,22,1,1,1,5])
C._addBC2Zone(a, 'overlap','BCOverlap',      [22,30,1,1,1,5])
C._addBC2Zone(a, 'wallOut','BCWallViscous',  [30,200,1,1,1,5])
t_curvi = C.newPyTree(['Base', a])

zones = Internal.getZones(t_curvi)
#les zones curviligne possedanr raccord chimere avec zone Cart IBC doit avoir la racine "joinIBC" dans leur nom 
zones[0][0]='curvi_joinIBC'
t_curvi = X.connectMatchPeriodic(t_curvi, translation=[0.,0.,0.04])
#stretch maillage plaque direction normal paroi
for z in zones:
  coordy =  Internal.getNodeFromName(z,'CoordinateY')[1]
  sh = np.shape(coordy)
  for j in range(1,sh[1]):
     coordy[:,j,:]=coordy[:,j-1,:]+0.002*1.02**j


U0= 120.0
P0= 101325.0
L0= 1.
R0= 1.18
T0= 279.15

equation = 'NSLaminar'
#equation = 'NSTurbulent'
C._addState(t_curvi, 'GoverningEquations', equation )
C._addState(t_curvi, 'EquationDimension', 3)
C._addState(t_curvi, UInf=U0, RoInf=R0, PInf=P0, LInf=L0, alphaZ=0., adim='dim3')
C._initVars(t_curvi,"{centers:Density}=1.18")
C._initVars(t_curvi,"{centers:VelocityX}=0.01")
C._initVars(t_curvi,"{centers:VelocityY}=0")
C._initVars(t_curvi,"{centers:VelocityZ}=0")
C._initVars(t_curvi,"{centers:Temperature}=279.15")

C._initVars(t_3d,"{centers:Density}=1.18")
C._initVars(t_3d,"{centers:VelocityX}=0.01")
C._initVars(t_3d,"{centers:VelocityY}=0")
C._initVars(t_3d,"{centers:VelocityZ}=0")
C._initVars(t_3d,"{centers:Temperature}=279.15")

#Calcul interpolation entre maillage curviligne  et cartesien IBC et ajout Ghost  au maillage curvi
t_final, tc_final=  App.setInterpData_Hybride(t_3d, tc_3d, t_curvi)

numb={"temporal_scheme": "implicit", "ss_iteration":10, "omp_mode":1,"modulo_verif":25}
numz={"time_step": 0.0001, "scheme":"roe","slope":"minmod", "time_step_nature":"global", "cfl":1.}

Fast._setNum2Zones(t_final, numz); Fast._setNum2Base(t_final, numb)

#tc_final=None
(t_final, tc_final, metrics) = FastS.warmup(t_final, tc_final,verbose=1)


nit = 100; time = 0.
for it in range(nit):
    FastS._compute(t_final, metrics, it, tc_final)
    if it%25==0:
       FastS.displayTemporalCriteria(t_final, metrics, it)

#C.convertPyTree2File(t_final, 'restart.cgns')
Internal._rmNodesByName(t_final, '.Solver#Param')
Internal._rmNodesByName(t_final, '.Solver#ownData')
Internal._rmNodesByName(t_final, '.Solver#dtloc')
Internal._rmNodesByName(t_final, '*P1*')
Internal._rmNodesByName(t_final, '*M1*')
Internal._rmNodesFromType(t_final, 'Rind_t')
test.testT(t_final, 2)

