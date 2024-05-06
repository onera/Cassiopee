# - moveCamera (pyTree) -
import Geom.PyTree as D
import CPlot.PyTree as CPlot
import Converter.PyTree as C
import Transform.PyTree as T
import Generator.PyTree as G

# Model
a = D.sphere((0,0,0), 1., N=20)
a = C.convertArray2Hexa(a); a = G.close(a)
CPlot.display(a, posCam=(3,-1,0.7), posEye=(0,0,0))

t = 0.
for i in range(1000):
    # change model
    C._initVars(a, '{df}=0.1*cos(%f)*sin(10*pi*{CoordinateX})'%(t))
    b = T.deformNormals(a, 'df')
    # Move camera
    CPlot.moveCamera([(3,-1,0.7),(3,5,0.7),(3,7,0.7)], N=1000, pos=i)
    CPlot.display(b)
    t += 0.05
