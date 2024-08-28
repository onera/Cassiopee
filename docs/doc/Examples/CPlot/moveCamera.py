# - moveCamera (array) -
import Geom as D
import Converter as C
import Transform as T
import CPlot

# Model
a = D.sphere((0,0,0), 1., N=20)
CPlot.display(a, posCam=(3,-1,0.7), posEye=(0,0,0))

t = 0.
for i in range(1000):
    # change model
    defo = C.initVars(a, '{df}=0.1*cos(%f)*sin(10*pi*{x})'%(t))
    defo = C.extractVars(defo, ['df'])
    b = T.deformNormals(a, defo)
    # Move camera
    CPlot.moveCamera([(3,-1,0.7),(3,5,0.7),(3,7,0.7)], N=1000, pos=i)
    CPlot.display(b)
    t += 0.05
