# - delete (array) -
import Generator as G
import CPlot
import time

a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
b = G.cart( (11,0,0), (1,1,1), (10,10,10) )

CPlot.display([a,b]); time.sleep(1)
CPlot.delete([0]); CPlot.render(); time.sleep(1)
CPlot.delete(['S-Zone 1']); CPlot.render(); time.sleep(1)
