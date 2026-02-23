# - display (array) -
import Generator as G
import CPlot
import Converter as C

def F(x,y): return x*x + y*y

a = G.cart((0,0,0),(1,1,1),(18,28,1))
a = C.initVars(a, 'F', F, ['x','y'])

CPlot.display(a, mode=3, scalarField=0, dim=2,
              displayIsoLegend=1)
