# - loadView (pyTree) -
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import Generator.PyTree as G

a = G.cart( (0,0,0), (1,1,1), (10,10,1) )
t = C.newPyTree(['Base', 2, a])

# isoScales specify field no, niso for this field, min, max for this field
t = CPlot.addRender2PyTree(t, slot=0, posCam=(1,1,1), posEye=(10,1,1),
                           mode='Solid', niso=10, isoScales=[0, 10, 0., 1.],
                           isoEdges=0.1, isoLight=1, colormap='Blue2Red')

CPlot.display(t)
import time; time.sleep(0.1)
CPlot.loadView(t, slot=0)
