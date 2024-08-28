# - getActivePointIndex (array) -
import Generator as G
import CPlot
import time

a = G.cartTetra( (0,0,0), (1,1,1), (5,5,1) )
CPlot.display([a], dim=2)

l = []
while l == []:
    l = CPlot.getActivePointIndex(); time.sleep(0.1)
print('ActivePointIndex : ', l)
#>> ActivePointIndex :  [16, 19, 3, 0, 0]
