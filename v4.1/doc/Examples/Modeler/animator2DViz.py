# - Animator2D -
# Permet de visualiser les differents billboards dans une image
# z: zoom, s: dezoom, arrows: move
import Modeler.Animator2D as Animator2D
import time
import CPlot

# Fichier image
FILE = 'Images/monstres.png'
# Nbre de colonnes
Ni = 1
# Nbre de lignes
Nj = 3

a = Animator2D.Animator2D()

# Always register images before opening display
a.registerImage('pix', FILE, Ni, Nj, speed=1./(Ni*Nj))

# Open display
a.openDisplay()

# Draw images
for j in range(Nj):
    for i in range(Ni):
        name = '%d-%d'%(i,j)
        a.drawImage('pix'+name, imageKey='pix', pos=(0.3*i-0.3,0.3*j-0.3,0), scale=0.3)
        a.setAnimation('pix'+name, i,j, imageKey='pix')
        a.drawText('txt'+name, pos=(0.3*i-0.3,0.3*j-0.3,0.01), text=name, scale=0.002)

posCam = CPlot.getState('posCam')
posEye = CPlot.getState('posEye')

dep = 0.05

# Move
while 1 != 2:
    l = ''
    while l == '':
        time.sleep(0.01); l = CPlot.getKeyboard();
    CPlot.resetKeyboard()
    v = ord(l[0])
    if v == 3: # droite
        posCam = (posCam[0]-dep,posCam[1],posCam[2])
        posEye = (posEye[0]-dep,posEye[1],posEye[2])
        CPlot.setState(posCam=posCam, posEye=posEye)
    elif v == 4:
        posCam = (posCam[0]+dep,posCam[1],posCam[2])
        posEye = (posEye[0]+dep,posEye[1],posEye[2])
        CPlot.setState(posCam=posCam, posEye=posEye)
    elif v == 1:
        posCam = (posCam[0],posCam[1]+dep,posCam[2])
        posEye = (posEye[0],posEye[1]+dep,posEye[2])
        CPlot.setState(posCam=posCam, posEye=posEye)
    elif v == 2:
        posCam = (posCam[0],posCam[1]-dep,posCam[2])
        posEye = (posEye[0],posEye[1]-dep,posEye[2])
        CPlot.setState(posCam=posCam, posEye=posEye)
    elif l[0] == 'z':
        delta = abs(posCam[2]-posEye[2])
        if delta > 0.01:
            dep = delta*0.05
            posCam = (posCam[0],posCam[1],posCam[2]-delta*0.1)
            posEye = (posEye[0],posEye[1],posEye[2])
            CPlot.setState(posCam=posCam, posEye=posEye)
    elif l[0] == 's':
        if delta < 10:
            delta = abs(posCam[2]-posEye[2])
            dep = delta*0.05
            posCam = (posCam[0],posCam[1],posCam[2]+delta*0.1)
            posEye = (posEye[0],posEye[1],posEye[2])
            CPlot.setState(posCam=posCam, posEye=posEye)
