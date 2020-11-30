# - evalPosition (pyTree) - 
import RigidMotion.PyTree as R
import KCore.test as test

def F(a):
    a = R.setPrescribedMotion2(a, 'rotor', transl_speed=(0.1,0,0))
    b = R.evalPosition(a, time=0.1)
    return b

test.stdTestT(F)
