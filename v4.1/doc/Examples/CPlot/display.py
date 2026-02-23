# - display (array) -
import Generator as G
import CPlot
import Transform as T

a = G.cart((0,0,0),(1,1,1),(18,28,3))
CPlot.display(a, mode='mesh')

for i in range(360):
    a = T.rotate(a, (9, 14, 3.5), (0,0,1), 1.)
    CPlot.display(a)
