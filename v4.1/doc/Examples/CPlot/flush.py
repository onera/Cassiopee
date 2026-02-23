# - flush (array) -
import Generator as G
import CPlot
import Transform as T
import time

a = G.cart((0,0,0),(1,1,1),(180,280,300))
CPlot.display(a, mode=0)
time.sleep(2)
#CPlot.flush()
time.sleep(2)
a = T.rotate(a, (0,0,0), (0,0,1), 20.)
CPlot.display(a, mode=0)
CPlot.cplot.show()
