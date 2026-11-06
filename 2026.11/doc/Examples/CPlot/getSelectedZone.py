# - getSelectedZone (array) -
import Generator as G
import CPlot
import time

a = G.cart((0,0,0), (1,1,1), (5,5,5))
CPlot.display(a)

nz = -1
while nz == -1:
    nz = CPlot.getSelectedZone(); time.sleep(0.1)
print('One zone has been selected: %d.'%nz)
