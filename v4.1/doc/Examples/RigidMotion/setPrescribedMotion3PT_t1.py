# - setPrescribedMotion3 (pyTree) -
# Motion defined by a constant rotation and translation speed
import RigidMotion.PyTree as R
import Geom.PyTree as D
import KCore.test as test

a = D.sphere((1.2,0.,0.), 0.2, 30)
a = R.setPrescribedMotion3(a, 'mot', transl_speed=(1,0,0),
                           axis_pnt=(1.,2.,0.), axis_vct=(0.,0.,1.),
                           omega=0.2)
test.testT(a, 1)
