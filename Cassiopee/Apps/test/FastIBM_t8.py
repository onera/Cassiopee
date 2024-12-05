# - Fast.IBM -
# NS, para, frontType=2
import Apps.Fast.IBM as App
import Converter.PyTree as C
import Converter.Internal as Internal
import math
import numpy as np
import KCore.test as test
test.TOLERANCE = 1.e-6

LOCAL = test.getLocal()

myApp = App.IBM(format='single')
myApp.set(numb={"temporal_scheme": "implicit",
                "ss_iteration":3,
                "omp_mode":1})
myApp.set(numz={"time_step": 0.0007,
                "scheme":"roe_min",
                "time_step_nature":"local",
                "cfl":4.})

tb_2d = C.convertFile2PyTree('naca1DNS.cgns')

# Prepare
t_2d, tc_2d = App.prepare1('naca1DNS.cgns', t_out=LOCAL+'/t.cgns', tc_out=LOCAL+'/tc.cgns', frontType=2, cleanCellN=False)

Internal._rmNodesFromType(tc_2d,'Rind_t')
Internal._rmNodesFromName(tc_2d,Internal.__GridCoordinates__)
test.testT(tc_2d, 1)

#on renomme les zones
i=0
for z in Internal.getZones(t_2d):
    z[0]= "Cart."+str(i)+"X0"
    i += 1

## determine dx=dy for each zone & store per zone
dict_ZEXT={}
hmin = 1.e30
for z in Internal.getZones(t_2d):
    h = abs(C.getValue(z,'CoordinateX',0)-C.getValue(z,'CoordinateX',1))
    print("dx=",h)
    dict_ZEXT[z[0]]=h
    if h < hmin : hmin = h

## go from dx to dx/dx_min
Nlevels=1
for i in dict_ZEXT:
    dict_ZEXT[i]= math.log( int(dict_ZEXT[i]/hmin + 0.00000001)  , 2)
    if dict_ZEXT[i] +1  > Nlevels : Nlevels = dict_ZEXT[i] +1

## get number of levels
print("Nlevel, hmin", Nlevels, hmin)


## Create the dict with Nz for each zone
dictNz      = {}

Nz_min = 8
Nz_max = 8
NzLoc=np.empty(int(Nlevels), np.int32)

for l in range( int(Nlevels) ):
    NzLoc[l] = Nz_max
Nlevels_tg = math.log( Nz_max/Nz_min, 2 ) +1 

for z in Internal.getZones(t_2d):

    level = int( dict_ZEXT[ z[0] ] )
    print("Nz local",  z[0], level)
    dictNz[z[0]]=  NzLoc[level]
    print("Nz local", NzLoc[level], z[0], level)

extrusion ='cart'; span = 0.078
#extrusion ='cyl'; span = 22.5

t_3d, tb_3d = App.extrudeCartesian(t_2d, tb_2d, extrusion=extrusion, NPas=5, span=span, Ntranche=1, dictNz=dictNz)

interpDataType = 1 # on suppose maillage non cartesion pour interp
order          = 2

t_3d, tc_3d = App.prepare1(tb_3d, None, None, t_in=t_3d, extrusion=extrusion, interpDataType=interpDataType, order=order, vmin=11)

C.convertPyTree2File(t_3d, LOCAL+'/t.cgns')
C.convertPyTree2File(tc_3d, LOCAL+'/tc.cgns')

# Compute
t, tc = myApp.compute(LOCAL+'/t.cgns', LOCAL+'/tc.cgns', t_out=LOCAL+'/restart.cgns', tc_out=LOCAL+'/tc_restart.cgns', nit=100)
t = C.convertFile2PyTree(LOCAL+'/restart.cgns')
Internal._rmNodesByName(t, '.Solver#Param')
Internal._rmNodesByName(t, '.Solver#ownData')
Internal._rmNodesByName(t, '.Solver#dtloc')
Internal._rmNodesFromType(t, 'Rind_t')
test.testT(t, 2)

# Post
t, zw = myApp.post(tb_3d, LOCAL+'/restart.cgns', LOCAL+'/tc_restart.cgns', t_out=LOCAL+'/out.cgns', wall_out=LOCAL+'/wall.cgns')
C._rmVars(t, 'mutsmu') # pas assez precis
Internal._rmNodesFromType(t, 'Rind_t')
Internal._rmNodesByName(t, '.Solver#dtloc')
test.testT(t, 3)
