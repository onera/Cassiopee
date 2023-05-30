# - display (array) -
# display offscreen using OpenGL
import CPlot
import Transform as T
import Geom as D
import KCore.config
import KCore.test as test
import os

if KCore.config.CPlotOffScreen: 
    print("CPlot: no offscreen GL installed.")
    os._exit(1)

a = D.sphere((0,0,0), 1.0, N=200)
b = D.sphere((0,1,0), 0.8, N=200)

# direct export
CPlot.display(a, mode=0, offscreen=2, export='out.png')
CPlot.finalizeExport(2)
test.testF('out.png', 1)

# composite export
CPlot.display(a, mode=0, offscreen=3, export='out2.png')
CPlot.finalizeExport(3)
CPlot.display(b, mode=0, offscreen=4, export='out2.png')
CPlot.finalizeExport(4)
test.testF('out2.png', 2)

os._exit(0)
