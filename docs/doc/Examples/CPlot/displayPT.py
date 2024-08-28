# - display (pyTree) -
import Generator.PyTree as G
import CPlot.PyTree as CPlot
import Transform.PyTree as T

a = G.cart((0,0,0),(1,1,1),(18,28,3))

for i in range(360):
    a = T.rotate(a, (9, 14, 3.5), (0,0,1), 1.)
    CPlot.display(a)