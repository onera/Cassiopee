# - Fast.MB -
import Apps.Fast.MB as Apps_MB
import Apps.Fast.Common as Apps_Common
import KCore.test as test
import Converter.PyTree as C
import Converter.Internal as Internal
import Converter.Mpi as Cmpi
import FastC.PyTree as FastC
test.TOLERANCE = 2.e-5

LOCAL = test.getLocal()

if Cmpi.rank == 0: # prep en seq pour l'instant
    t, tc = Apps_MB.prepare('naca.cgns', t_out=LOCAL+'/t.cgns', tc_out=LOCAL+'/tc.cgns', NP=Cmpi.size)
    test.testT(tc, 2)
Cmpi.barrier()

if Cmpi.size > 1:
    import FastS.Mpi as FastS
    rank = Cmpi.rank; size = Cmpi.size
else:
    import FastS.PyTree as FastS
    rank = 0; size = 1

t,tc,ts,graph = FastC.load(LOCAL+'/t.cgns', LOCAL+'/tc.cgns', split='single')

# Numerics
numb={"temporal_scheme": "implicit",
      "ss_iteration":3}
numz={"time_step": 0.0007,
      "scheme":"roe_min",
      "time_step_nature":"local",
      "cfl":4.}

FastC._setNum2Zones(t, numz); FastC._setNum2Base(t, numb)

(t, tc, metrics) = FastS.warmup(t, tc, graph)

it0 = 0; time0 = 0.
first = Internal.getNodeFromName1(t, 'Iteration')
if first is not None: it0 = Internal.getValue(first)
first = Internal.getNodeFromName1(t, 'Time')
if first is not None: time0 = Internal.getValue(first)
time_step = Internal.getNodeFromName(t, 'time_step')
time_step = Internal.getValue(time_step)

if 'modulo_verif' in numb: moduloVerif = numb['modulo_verif']
else: moduloVerif = 200

NIT = 300
for it in range(NIT):
    FastS._compute(t, metrics, it, tc, graph)
    if it%moduloVerif == 0:
        if rank == 0: print('- %d / %d - %f'%(it+it0, NIT+it0, time0))
        FastS.display_temporal_criteria(t, metrics, it, format='double')
        #if it%50 == 0:
        #    import CPlot.PyTree as CPlot
        #    CPlot.display(t, dim=2, mode='Scalar', scalarField='Density')
    time0 += time_step

# time stamp
Internal.createUniqueChild(t, 'Iteration', 'DataArray_t', value=it0+NIT)
Internal.createUniqueChild(t, 'Time', 'DataArray_t', value=time0)

FastC.save(t, LOCAL+'/restart.cgns', split='single', compress=0)
if Cmpi.size > 1: Cmpi.barrier()

Cmpi.barrier()
if Cmpi.rank == 0:
    t = C.convertFile2PyTree(LOCAL+'/restart.cgns')
    Internal._rmNodesByName(t, '.Solver#Param')
    Internal._rmNodesByName(t, '.Solver#ownData')
    Internal._rmNodesByName(t, '.Solver#dtloc')
    test.testT(t, 1)
