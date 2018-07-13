# - moveCamera (array) -
import Geom as D
import CPlot

# Model
a = D.sphere((0,0,0), 1., N=20)
CPlot.display(a, posCam=(3,-1,0.7), posEye=(0,0,0))

# Move camera
CPlot.moveCamera([(3,-1,0.7),(3,5,0.7),(3,7,0.7)], N=1000, speed=1)
