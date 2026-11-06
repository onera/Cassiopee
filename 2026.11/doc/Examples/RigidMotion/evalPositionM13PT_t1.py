# - evalPositionM1 pour motion 3 (pyTree) -
# Motion defined by a constant rotation and translation speed
import RigidMotion.PyTree as R
import Geom.PyTree as D
import Converter.Internal as Internal
import KCore.test as test

a = D.sphere((1.2,0.,0.), 0.2, 30)
a = R.setPrescribedMotion3(a, 'mot', transl_speed=(0.1,0,0),
                           axis_pnt=(1.,2.,0.), axis_vct=(0.,0.,1.),
                           omega=0.2)

b = R.evalPosition(a, time=1.); b[0]='moved'
coords = [Internal.getNodeFromName(b, 'CoordinateX')[1],
          Internal.getNodeFromName(b, 'CoordinateY')[1],
          Internal.getNodeFromName(b, 'CoordinateZ')[1]]
coords0 = R.evalPositionM1(coords, b, time=1.)
c = Internal.copyRef(b); c[0] = 'back'
Internal.getNodeFromName(c, 'CoordinateX')[1] = coords0[0]
Internal.getNodeFromName(c, 'CoordinateY')[1] = coords0[1]
Internal.getNodeFromName(c, 'CoordinateZ')[1] = coords0[2]
test.testT([a,b,c], 1)
