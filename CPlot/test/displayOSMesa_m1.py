# - display (array) -
# display offscreen using OSMesa parallel
import CPlot
import Transform as T
import Geom as D
import KCore.config
import KCore.test as test

if not KCore.config.CPlotOffScreen: 
    print("CPlot: no offscreen osmesa installed.")
    import sys; sys.exit()

a = D.sphere((0,Cmpi.rank,0), 1, N=200)

# direct export
CPlot.display(a, mode=0, offscreen=7, export='out.png')
test.testF('out.png', 1)
