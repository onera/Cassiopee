# - Animator2D -
# Test selector
import Modeler.Animator2D as Animator2D

a = Animator2D.Animator2D()

# Open display
a.openDisplay()

a.drawSelector('selector', pos=(-0.1,0,0), w=0.2, h=0.1, e=0.01)
