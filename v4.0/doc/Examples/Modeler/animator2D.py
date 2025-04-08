# - Animator2D -
# Animation billboards d'un personnage

import Modeler.Animator2D as Animator2D
import time
import random

state = 0 # sourire

def cligne():
    a.erase('yeux1', render=False)
    a.drawImage('yeux2', pos=(0,0,0.01), scale=0.6)
    time.sleep(0.2)
    a.erase('yeux2', render=False)
    a.drawImage('yeux1', pos=(0,0,0.01), scale=0.6)

def souris():
    global state
    if state == 1:
        a.erase('bouche2', render=False)
        a.drawImage('bouche1', pos=(0,0,0.01), scale=0.6)
        state = 0

def triste():
    global state
    if state == 0:
        a.erase('bouche1', render=False)
        a.drawImage('bouche2', pos=(-0.005,0,0.01), scale=0.6)
        state = 1


a = Animator2D.Animator2D()

# Always register images before opening display
a.registerImage('perso', 'Images/personnage.png')
a.registerImage('bouche1', 'Images/bouche_souriante.png')
a.registerImage('bouche2', 'Images/bouche_non_souriante.png')
a.registerImage('yeux1', 'Images/yeux_ouverts.png')
a.registerImage('yeux2', 'Images/yeux_fermes.png')

# Open display
a.openDisplay()

# Draw images
a.drawImage('perso', pos=(0,0,0), scale=0.6)
a.drawImage('yeux1', pos=(0,0,0.01), scale=0.6)
a.drawImage('bouche1', pos=(0,0,0.01), scale=0.6)

time.sleep(2.)

# Animate
while 1 != 2:
    r = random.randint(0,20)
    if r == 0 or r == 1 or r == 2: cligne()
    elif r == 5:
        triste()
    elif r == 8: souris()
    time.sleep(0.5)
