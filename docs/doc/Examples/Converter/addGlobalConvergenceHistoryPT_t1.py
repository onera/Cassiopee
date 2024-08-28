# - addGlobalConvergenceHistory (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.elsAProfile as elsAProfile
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (10,10,10))
t = C.newPyTree(['Base',a])
tp = elsAProfile.addGlobalConvergenceHistory(t)
test.testT(tp, 1)
# update node with another norm value
tp = elsAProfile.addGlobalConvergenceHistory(tp,normValue=1)
test.testT(tp,2)
