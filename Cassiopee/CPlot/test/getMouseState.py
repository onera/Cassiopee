# - getMouseState (array) -
import Generator as G
import CPlot
import time

a = G.cartTetra( (0,0,0), (1,1,1), (5,5,1) )
CPlot.display([a], dim=2)

c = 1000
while c > 0:
    l = CPlot.getMouseState(); time.sleep(0.5)
    print(l)
    c -= 1
