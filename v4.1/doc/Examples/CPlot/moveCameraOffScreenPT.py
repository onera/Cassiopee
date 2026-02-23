# - moveCamera (pyTree) -
import Geom.PyTree as D
import CPlot.PyTree as CPlot
import Converter.PyTree as C
import Transform.PyTree as T
import Generator.PyTree as G

# Model
a = D.sphere((0,0,0), 1., N=20)
a = C.convertArray2Hexa(a); a = G.close(a)

t = 0.
for i in range(100):
    # change model
    C._initVars(a, '{df}=0.1*cos(%f)*sin(10*pi*{CoordinateX})'%(t))
    b = T.deformNormals(a, 'df')
    # Move camera
    (posCam, posEye, dirCam) = CPlot.moveCamera([(3,-1,0.7),(3,5,0.7),(3,7,0.7)], N=100, pos=i)
    CPlot.display(b, posCam=posCam, posEye=posEye, dirCam=dirCam, offscreen=1, export='image%i.png'%i)
    t += 0.05
