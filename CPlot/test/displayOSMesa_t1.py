# - display (array) -
# display offscreen using OSMesa
import CPlot
import Transform as T
import Geom as D
import KCore.config
import KCore.test as test

if not KCore.config.CPlotOffScreen: 
    print("CPlot: no offscreen osmesa installed.")
    import sys; sys.exit()

a = D.sphere((0,0,0), 1, N=200)
b = D.sphere((0,1,0), 0.8, N=200)

# direct export
CPlot.display(a, mode=0, offscreen=1, export='out.png')
CPlot.finalizeExport(1)
test.testF('out.png', 1)

# composite export
CPlot.display(b, mode=0, offscreen=5, export='out2.png')
CPlot.finalizeExport(5)
CPlot.display(b, mode=0, offscreen=6, export='out2.png')
CPlot.finalizeExport(6)
test.testF('out2.png', 2)
