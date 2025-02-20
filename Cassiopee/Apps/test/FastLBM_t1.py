# - PSEUDO ISENTROPIC VORTEX MULTIBLOCK (pyTree) -
import Apps.Fast.LBM as Apps_LBM
import Connector.PyTree as X
import Converter.Internal as Internal
import Converter.PyTree as C
import FastC.PyTree as FastC
import FastLBM.PyTree as FastLBM
import Generator.PyTree as G
import KCore.test as test
import Transform.PyTree as T
import math
import numpy

VARSMACRO = ['Density','VelocityX','VelocityY','VelocityZ','Temperature']

# Geometry and mesh
nit   = 100
NG    = 1
x_min = -0.5
dx    = 0.01
L     = 1.
Nx = 101;
Ny = Nx;
Nz = 9
c0=343.2
n = Nx-1

tree=dict()
tree['mach']=0.1
tree['maille']=0.01
tree['timestep']=tree['maille']/(c0*numpy.sqrt(3.))
tree['reynolds']=100.
tree['char_length']=L/10. # R0
tree['rho0']=1.

# Fluid
rho0 = 1.
u0 = tree["mach"]* c0

# Vortex
R0 = tree['char_length']
kappa = 0.14
dt = dx / (3. ** 0.5 * c0)
#-------------------------------------------------------
#-------------------------------------------------------
z_min = -4*dx
a1 = G.cart((x_min,x_min,z_min), (dx,dx,dx), (Nx,Ny,Nz))
x_0 = C.getMeanValue(a1,'CoordinateX')
y_0 = C.getMeanValue(a1,'CoordinateY')
z_0 = C.getMeanValue(a1,'CoordinateZ')
zmean = C.getMeanValue(a1,'CoordinateZ')
a1 = T.splitNParts(a1,8)
t = C.newPyTree(['Base',a1])
t,tc = Apps_LBM.prepare(t, t_out=None, tc_out=None, NP=0, translation=[(Nx-1)*dx, (Ny-1)*dx,(Nz-1)*dx],NG=NG)

#-------------------------
# Initialization
#-------------------------
rho0Kap2 = rho0*kappa**2
kapC0=kappa*c0
R0_2Inv= 1./(R0*R0)
v_adim=1./(math.sqrt(3)*c0)

R0 = L/10.; R0_2Inv= 1./(R0*R0); R0_3Inv= 1./(R0*R0*R0)
epsilon = 0.07*c0; eps_2 = epsilon**2
Tref = 1.0; Rhoref = 1.0

C._initVars(t,"{centers:r2}=({centers:CoordinateX}-%g)**2+({centers:CoordinateY}-%g)**2"%(x_0,y_0))
C._initVars(t,'{centers:Density}= %20.16g*exp(-%20.16g/(2.*%20.16g**2)*exp(-{centers:r2}*%20.16g))'%(Rhoref,eps_2,c0,R0_2Inv))
C._initVars(t,'{centers:VelocityX}=%20.16g-%20.16g*(({centers:CoordinateY}-%20.16g)*%20.16g)*exp(-0.5*{centers:r2}*%20.16g)'%(u0,epsilon,y_0,1./R0,R0_2Inv))
C._initVars(t,'{centers:VelocityY}=%20.16g+%20.16g*(({centers:CoordinateX}-%20.16g)*%20.16g)*exp(-0.5*{centers:r2}*%20.16g)'%(0.,epsilon,x_0,1./R0,R0_2Inv))
C._initVars(t,'{centers:VelocityZ}=0.')
C._initVars(t,'{centers:Temperature}=%20.16g'%(Tref))

#-------------------------
# Adimensionnement
#-------------------------
C._initVars(t,'{centers:VelocityX}={centers:VelocityX}*%g'%v_adim)
C._initVars(t,'{centers:VelocityY}={centers:VelocityY}*%g'%v_adim)
C._initVars(t,"{centers:TurbulentDistance} =1.")
C._initVars(t,"{centers:cellN_IBC_LBM_1} =0.")
C._initVars(t,"{centers:cellN_IBC_LBM_2} =0.")
C._initVars(t,"{centers:cellN_IBC_LBM_3} =0.")
for v in VARSMACRO: C._cpVars(t,'centers:'+v,tc,v)
X._setInterpTransfers(t,tc,variables=VARSMACRO)
#-------------------------
# Compute
#-------------------------
# Numerics
numb                    = {'temporal_scheme':'explicit', 'ss_iteration':20,'omp_mode':0}
numz                    = {'scheme':'ausmpred'}
numz['cache_blocking_I']=1000000
numz['cache_blocking_J']=1000000
numz['cache_blocking_K']=1000000
numz['cache_blocking_I']=1000000

colop_select     = 'RR'
NQ_local         = 'D3Q19'
nu_local         = 0
numz["time_step"]           = 1
numz["LBM_velocity_set"]    = NQ_local
numz["lbm_c0"]              = math.sqrt(1./3.)
numz["lbm_dif_coef"]        = nu_local

numz["LBM_coll_model"]      =colop_select
numz["LBM_relax_time"]      =0.5


FastC._setNum2Zones(t, numz); FastC._setNum2Base(t, numb)

(t, tc, metrics)  = FastLBM.warmup(t, tc,nghost=NG,flag_initprecise=0)

for it in range(1,nit+1):
    if it%10==0:
        print("--------- iteration %d -------------"%it)
    FastLBM._compute(t, metrics, it, tc,layer='Python',nittotal=nit)

Internal._rmNodesByName(t, '.Solver#Param')
Internal._rmNodesByName(t, '.Solver#ownData')
# POST
Internal._rmNodesByName(t,'*M1*')
Internal._rmNodesByName(t,'*_P1')
Internal._rmNodesByName(t,'*cell*')
Internal._rmNodesByName(t,'*Qstar*')
Internal._rmNodesByName(t,'*Qeq*')
Internal._rmNodesByName(t,'*Qneq*')
Internal._rmNodesByName(t,'*SpongeCoef')

# reconstruction des variables macro
v_adiminv = math.sqrt(3)*c0
C._initVars(t,"{centers:VelocityX}={centers:VelocityX}*%g"%v_adiminv)
C._initVars(t,'{centers:VelocityY}={centers:VelocityY}*%g'%v_adiminv)

#C.convertPyTree2File(t,"restart.cgns")
test.testT(t,1)
