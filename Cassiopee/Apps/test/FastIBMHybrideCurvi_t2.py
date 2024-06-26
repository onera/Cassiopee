# - Fast.IBM -
import Apps.Fast.IBM as App
import Connector.PyTree as X
import Connector.IBM as X_IBM
import Converter.Internal as Internal
import Converter.PyTree as C
import Fast.PyTree as Fast
import FastS.PyTree as FastS
import Generator.PyTree as G
import Geom.IBM    as D_IBM
import Geom.PyTree as D
import Initiator.PyTree as Initiator
import KCore.test as test
import math
import numpy as np
test.TOLERANCE = 1.e-6

LOCAL = test.getLocal()

#Fabrication case 2d pour maillage octree 
lines=[]
lines.append( D.line((0.05,-0.15, 0), (0.1 ,-0.08,0), N=10) )
lines.append( D.line((0.10,-0.03, 0), (0.1 ,-0.08,0), N=10) )
lines.append( D.line((0.1 ,-0.03, 0), (0.1 , 0.  ,0), N=8) )

lines.append( D.line((0.15,-0.08, 0), (0.15,-0.03,0), N=10) )
lines.append( D.line((0.15,-0.08, 0), (0.2 ,-0.15,0), N=10) )
lines.append( D.line((0.05,-0.15, 0), (0.2 ,-0.15,0), N=10) )

lines.append( D.line((0.15,-0.03, 0), (0.15, 0   ,0), N=8 ) )
lines.append( D.line((0.10, 0.01, 0), (0.15, 0.01,0), N=10) )
lines.append( D.line((0.15, 0.01, 0), (0.15, 0   ,0), N=3) )
lines.append( D.line((0.10, 0.01, 0), (0.10, 0   ,0), N=3) )
case = C.newPyTree(['BOBY', lines])

D_IBM._setSnear(case, 0.002)
D_IBM._setIBCType(case, "Musker")
D_IBM._setDfar(case, 0.)
D_IBM._setFluidInside(case)

zones = Internal.getZones(case)
for z in zones:
 if '8' in z[0] or '9'  in z[0]  or '6'  in z[0]:
    D_IBM._setIBCType(z, "overlap")
    D_IBM._setSnear(z, 0.008)
 if 'line.1' == z[0] or '10'  in z[0]:
    D_IBM._setSnear(z, 0.008)
 if '2' in z[0] or '0'  in z[0]:
    D_IBM._setSnear(z, 0.004)

zones[6][0]='line10'
zones[8][0]='line9'
U0= 120.0
P0= 101325.0
L0= 1.
R0= 1.18
T0= 279.15

equation = 'NSLaminar'
C._addState(case, 'GoverningEquations', equation )
C._addState(case, 'EquationDimension', 2)
C._addState(case, UInf=U0, RoInf=R0, PInf=P0, LInf=L0, alphaZ=0., adim='dim3')
#C.convertPyTree2File(case,'verifcase.cgns')

# Prepare cas 2d IBM
t_2d, tc_2d = App.prepare1(case, None, None, frontType=42 ,cleanCellN=False)

test.testT(tc_2d, 1)
test.testT(t_2d , 2)

#extrusion 3D IBM
extrusion ='cart'; span = 0.01

t_3d, tb_3d = App.extrudeCartesian(t_2d, case, extrusion=extrusion, NPas=5, span=span)

test.testT(tb_3d, 3)
test.testT(t_3d , 4)

#calcul interp IBM 3D; celln, distnce paroi 3D,.... deja obtenue lors de l'extrusion
interpDataType = 1 # on suppose maillage non cartesion pour interp
order          = 2
t_3d, tc_3d = App.prepare1(tb_3d,None, None, t_in=t_3d, extrusion=extrusion, interpDataType=interpDataType, order=order)

Internal._rmNodesFromType(tc_2d,'Rind_t')
Internal._rmNodesFromName(tc_2d,Internal.__GridCoordinates__)
test.testT(tc_3d, 5)
test.testT(t_3d , 6)


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
#t_final, tc_final=  App.setInterpData_Hybride(t_3d, tc_3d, t_curvi)
t_final, tc_final=  X_IBM.setInterpDataHybrid(t_3d, tc_3d, t_curvi)
test.testT(t_final , 7)
test.testT(tc_final, 8)

