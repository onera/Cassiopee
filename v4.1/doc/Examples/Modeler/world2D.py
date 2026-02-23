# - world2D -
import Modeler.World2D as World2D
import Modeler.Perso2D as Perso2D
import CPlot.PyTree as CPlot
import time

w = World2D.World('Images/level.txt')
w.render()
time.sleep(1.)

p = Perso2D.perso('hero', 'Images/serge2.png')
m = Perso2D.perso('serge', 'Images/serge.png')

p.build1(w)
p.place(w)
p.render1(w, True)

m.baseSpeed=50.
m.build1(w)
m.place(w, 2,2)
m.render1(w, True, False)

Perso2D.animator.prepare()
CPlot.render()

# Ready?
time.sleep(1.)

# Main loop
while 1 != 2:
    key = CPlot.getKeyboard()
    key = key[-1:]
    try: v = ord(key)
    except: v = key

    if v == 1: # up
        p.moveUp()
        CPlot.resetKeyboard()
        p.move(w)
        ret = p.checkCollision(w)
        p.render1(w); CPlot.render()
    elif v == 2: # down
        p.moveDown()
        CPlot.resetKeyboard()
        p.move(w)
        ret = p.checkCollision(w)
        p.render1(w); CPlot.render()
    elif v == 3: # left
        p.moveLeft()
        CPlot.resetKeyboard()
        p.move(w)
        ret = p.checkCollision(w)
        p.render1(w); CPlot.render()
    elif v == 4: # right
        p.moveRight()
        CPlot.resetKeyboard()
        p.move(w)
        ret = p.checkCollision(w)
        p.render1(w); CPlot.render()

    if v == 1 or v == 2 or v == 3 or v == 4:
        CPlot.resetKeyboard()
        p.move(w)
        ret = p.checkCollision(w)
        p.render1(w); CPlot.render()
        time.sleep(0.1)

    if v == 5 or v == 6 or v == 7 or v == 8 or v == '':
        # Finish animation
        cur = p.getCur()*16.
        if int(cur)/2-cur*0.5 == 0:
            if cur >= 1. and cur <= 4.: p.moveUp()
            if cur >= 5. and cur <= 8.: p.moveRight()
            if cur >= 9. and cur <= 12.: p.moveLeft()
            if cur >= 13. and cur <= 16.: p.moveDown()
            p.move(w)
            ret = p.checkCollision(w)
            p.render1(w); CPlot.render()
            time.sleep(0.2)
            CPlot.render()
    # NP mouvement
    if int(w.time*10) % 50000 == 0:
        m.autoMove()
        m.move(w)
        ret = m.checkCollision(w)
        m.render1(w, False, False)
        CPlot.render()

    w.advanceTime()
