# - travel* (array) -
import Geom as D
import Generator as G
import CPlot
import Converter as C
import Transform as T

# Model
a = D.sphere((0,0,0), 1., N=20)
a = C.convertArray2Hexa(a); a = G.close(a)
CPlot.display(a, posCam=(3,0,0), posEye=(0,0,0))

t = 0.
for i in range(1300):
    # change model
    defo = C.initVars(a, '{df}=0.1*cos(%f)*sin(10*pi*{x})'%(t))
    defo = C.extractVars(defo, ['df'])
    b = T.deformNormals(a, defo)
    CPlot.display(b)
    t += 0.05

    if i < 200: CPlot.travelLeft(0.001, N=3)
    elif i < 400: CPlot.travelRight(0.001, N=3)
    elif i < 600: CPlot.travelUp(0.001, N=3)
    elif i < 800: CPlot.travelIn(0.001, N=3)
