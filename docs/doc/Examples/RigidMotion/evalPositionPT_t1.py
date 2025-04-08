# - evalPosition (pyTree) -
import RigidMotion.PyTree as R
import KCore.test as test

def F(a):
    a = R.setPrescribedMotion1(a, 'trans', tx="{t}")
    b = R.evalPosition(a, time=0.1)
    return b

test.stdTestT(F)
