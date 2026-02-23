# - evalPosition (pyTree) -
import RigidMotion.PyTree as R
import KCore.test as test
def F(a):
    R._setPrescribedMotion3(a, 'constant', transl_speed=(0.1,0,0),
                            axis_pnt=(0.,0.,0.), axis_vct=(0,0,1), omega=1.)
    b = R.evalPosition(a, time=0.1)
    return b

test.stdTestT(F)
