# - setInterpData (pyTree) -
# case without cellN field: the whole receiver zone is interpolated
import Converter.PyTree as C
import Generator.PyTree as G
import Connector.PyTree as X
import Geom.PyTree as D
import KCore.test as test

# Tetra donor zone
a = G.cartTetra((0,0,-0.2),(0.01,0.01,0.1),(101,101,5))
pts = D.circle((0.5,0.5,0),0.05,N=20)
C._initVars(pts, 'cellN', 2); C._initVars(pts, 'centers:cellN',2)
# options to combine
notest = 1
for location in ['nodes', 'centers']:
    for stk in ['direct', 'indirect']:
        for pen in [0,1]:
            for nat in [0,1]:
                for order in [2,3,5]:
                    pts2 = X.setInterpData(pts, a, order=order, penalty=pen,
                                           nature=nat, loc=location,
                                           storage=stk, hook=None)
                    test.testT(pts2, notest)
                    notest += 1
