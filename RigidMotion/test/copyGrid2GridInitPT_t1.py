# - copyGrid2GridInit (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import RigidMotion.PyTree as R
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (10,10,10))
R._setPrescribedMotion3(a, 'motion', transl_speed=(1,0,0))
R._copyGrid2GridInit(a)
test.testT(a, 1)
