# - setSelectedZones (array) -
import Generator as G
import CPlot
import time

a1 = G.cart( (0,0,0), (1,1,1), (5,5,5) )
a2 = G.cart( (7,0,0), (1,1,1), (3,3,3) )
CPlot.display([a1, a2])

time.sleep(1.)
CPlot.setSelectedZones([(0,1), (1,1)])
