# - evalGridSpeed pour motion 3 (pyTree) -
# Motion defined by a constant rotation and translation speed
import RigidMotion.PyTree as R
import Converter.PyTree as C
import Geom.PyTree as D

a = D.sphere((1.2,0.,0.), 0.2, 30)
a = R.setPrescribedMotion3(a, 'mot', transl_speed=(0.2,0,0),
                           axis_pnt=(1.,2.,0.), axis_vct=(0.,0.,1.),
                           omega=0.2)

b = R.evalPosition(a, time=1.); b[0]='moved'
# must work on coordinates at time=1.
R._evalGridSpeed(b, time=1.)

# Comparing with finite difference
# must work on coordinates at time=0.
R._evalGridSpeed2(a, time=1.)
R._evalPosition(a, time=1.)

C.convertPyTree2File([a,b], 'out.cgns')
