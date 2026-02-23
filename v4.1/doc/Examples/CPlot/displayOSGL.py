# - display (array) -
# display offscreen using OpenGL
import CPlot
import Geom as D
import KCore.config
import KCore.test as test
import os

LOCAL = test.getLocal()

a = D.sphere((0,0,0), 1.0, N=200)
b = D.sphere((0,1,0), 0.8, N=200)

# direct export
CPlot.display(a, mode=0, offscreen=2, export=LOCAL+'/out.png')
CPlot.finalizeExport(2)
test.testF(LOCAL+'/out.png', 1)

# composite export
CPlot.display(a, mode=0, offscreen=3, export=LOCAL+'/out2.png')
CPlot.finalizeExport(3)
CPlot.display(b, mode=0, offscreen=4, export=LOCAL+'/out2.png')
CPlot.finalizeExport(4)
test.testF(LOCAL+'/out2.png', 2)

os._exit(0)