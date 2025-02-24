# - Fast.MB -
import Apps.Fast.MB as Apps_MB
import FastC.PyTree as FastC
import FastS.PyTree as FastS
import Converter.PyTree as C
import Converter.Internal as Internal
import KCore.test as test

LOCAL = test.getLocal()

t, tc = Apps_MB.prepare('naca.cgns', t_out=LOCAL+'/t.cgns', tc_out=LOCAL+'/tc.cgns', NP=0)
test.testT(tc, 2)

## Compute
numb={
    "temporal_scheme": "implicit",
    "ss_iteration":3
}
numz={
    "time_step": 0.0007,
    "scheme":"roe_min",
    "time_step_nature":"local",
    "cfl":4.
}

it0 = 0; time0 = 0.
FastC._setNum2Zones(t, numz); FastC._setNum2Base(t, numb)

(t, tc, metrics) = FastS.warmup(t, tc)

nit = 300
for it in range(nit):
    FastS._compute(t, metrics, it, tc)
    time0 += numz['time_step']

# time stamp
Internal.createUniqueChild(t, 'Iteration', 'DataArray_t', value=it0+nit)
Internal.createUniqueChild(t, 'Time', 'DataArray_t', value=time0)

FastC.save(t, LOCAL+'/restart.cgns', split='single', compress=0)

t = C.convertFile2PyTree(LOCAL+'/restart.cgns')
Internal._rmNodesByName(t, '.Solver#Param')
Internal._rmNodesByName(t, '.Solver#ownData')
test.testT(t, 1)
