# - getActivePoint (array) -
import Generator as G
import CPlot
import time

a = G.cart( (0,0,0), (1,1,1), (5,5,5) )
CPlot.display([a])

l = []
while l == []:
    l = CPlot.getActivePoint(); time.sleep(0.1)
print('ActivePoint: ', l)
#>> ActivePoint:  [3.9996489035268743, 2.127948294736359, 2.41771355073051]
