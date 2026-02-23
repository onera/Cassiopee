# - addPeriodicDataInSolverParam (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.elsAProfile as elsAProfile
import KCore.test as test
a = G.cylinder((0.,0.,0.), 0.1, 1., 0., 90., 5., (11,11,11))
t = C.newPyTree(['Base',a])
tp = elsAProfile.addPeriodicDataInSolverParam(t,rotationAngle=[0.,0.,90.],isChimera=True)
test.testT(tp,1)

a = G.cylinder((0.,0.,0.), 0.1, 1., 0., 90., 5., (11,11,11))
t = C.newPyTree(['Base',a])
tp = elsAProfile.addPeriodicDataInSolverParam(t,NAzimutalSectors=4,isChimera=True)
test.testT(tp,2)

a = G.cylinder((0.,0.,0.), 0.1, 1., 0., 90., 5., (11,11,11))
t = C.newPyTree(['Base',a])
elsAProfile._addPeriodicDataInSolverParam(t,NAzimutalSectors=4,isChimera=True)
test.testT(t,3)
