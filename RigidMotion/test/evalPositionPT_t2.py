# - evalPosition (pyTree) - 
import RigidMotion.PyTree as R
import KCore.test as test

# Utilise le RigidMotion du solveur
try:
    import Cassiopee as K
    import elsA_user as E
except: # no error in this case
    import sys
    sys.exit()

def F(a):
    a = R.setPrescribedMotion2(a, 'rotor', transl_speed=(0.1,0,0))
    b = R.evalPosition(a, time=0.1)
    return b

test.stdTestT(F)
