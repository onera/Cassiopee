# - delete (array) -
import Generator as G
import CPlot
import time

a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
b = G.cart( (0,10,0), (1,1,1), (10,10,10) )
c = G.cartTetra( (11,0,0), (1,1,1), (10,10,10) )

CPlot.display([a,b,c]); time.sleep(1)
CPlot.delete([1]); CPlot.render(); time.sleep(1)
