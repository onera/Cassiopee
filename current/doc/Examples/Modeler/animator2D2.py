# - Animator2D -
# Animation billboard avec plusieurs billboards dans une image
import Modeler.Animator2D as Animator2D
import time

a = Animator2D.Animator2D()

# Always register images before opening display
a.registerImage('serge', 'Images/serge.png', 4, 4, speed=1./16.)

# Open display
a.openDisplay()

# Draw images
a.drawImage('serge', pos=(0,0,0), scale=0.5)

time.sleep(2.)

# Animate
while 1 != 2:
    a.animate('serge')
    time.sleep(0.5)
