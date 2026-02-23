# - connectMatch (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.elsAProfile as elsAProfile

a = G.cylinder((0.,0.,0.), 0.1, 1., 0., 90., 5., (11,11,11))
t = C.newPyTree(['Base', a])
tp = elsAProfile.addPeriodicDataInSolverParam(t, rotationAngle=[0.,0.,90.], isChimera=True)
C.convertPyTree2File(tp, 'out.cgns')
