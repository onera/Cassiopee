# - display (array) -
# display offscreen using OSMesa
import CPlot
import Geom as D

a = D.sphere((0,0,0), 1, N=200)

CPlot.display(a, mode=0, offscreen=1, export='out.png')
