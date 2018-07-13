# - moveCamera (pyTree) -
import Geom.PyTree as D
import CPlot.PyTree as CPlot

# Model
a = D.sphere((0,0,0), 1., N=20)
CPlot.display(a)
CPlot.setState(posCam=(1.5,0,0), posEye=(0,0,0))

CPlot.moveCamera([(1.5,0,0),(0.75,0.75,0),(0,0,0)])
