# - Animator2D -
# Essai particles
import Modeler.Animator2D as Animator2D
import time

a = Animator2D.Animator2D()

# Open display
a.openDisplay()

a.drawParticles('particules1', pos=(0,0,0), Np=30)

for i in range(100):
    a.animateParticles('particules1', 0.001, 0.1)
    time.sleep(0.1)
