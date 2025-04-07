# - travel* (pyTree) -
import Geom.PyTree as D
import Generator.PyTree as G
import CPlot.PyTree as CPlot
import Converter.PyTree as C
import Transform.PyTree as T

# Model
a = D.sphere((0,0,0), 1., N=20)
a = C.convertArray2Hexa(a); a = G.close(a)
CPlot.display(a, posCam=(3,0,0), posEye=(0,0,0))

time = 0.
for i in range(1300):
    # change model
    C._initVars(a, '{df}=0.1*cos(%f)*sin(10*pi*{CoordinateX})'%(time))
    b = T.deformNormals(a, 'df')
    CPlot.display(b)
    time += 0.05

    if i < 200: CPlot.travelLeft(0.001, N=3)
    elif i < 400: CPlot.travelRight(0.001, N=3)
    elif i < 600: CPlot.travelUp(0.001, N=3)
    elif i < 800: CPlot.travelIn(0.001, N=3)
