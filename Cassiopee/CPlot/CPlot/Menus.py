# Create menus and submenus in CPlot
from . import PyTree as CPlot
import time
dt = 0.01

# Cree un submenu, affiche les items et retourne l'item selectionne
# IN: items: liste de chaines
def submenu(items):
    import Geom.PyTree  as D
    import Transform.PyTree as T
    import Converter.PyTree as C
    nitems = len(items)
    out =  []; c = 0
    for i in items:
        a = D.text2D(i)
        CPlot._addRender2Zone(a, color='White')
        T._translate(a, (0,-c,0))
        out.append(a)
        c += 10
    t = C.newPyTree(['Base']); t[2][1][2] += out
    CPlot.display(t, mode='Render')
    CPlot.setState(activateShortCuts=0)
    current = 0
    go = True; l = ''
    while go:
        while l == '':
            l = CPlot.getKeyboard(); time.sleep(dt)
        v = ord(l[0]); print(v)
        if v == 1: # up
            b = CPlot.addRender2Zone(out[current], color='White')
            CPlot.replace(t, 1, current, b)
            current -= 1
            if current < 0: current = nitems-1
            b = CPlot.addRender2Zone(out[current], color='Red')
            CPlot.replace(t, 1, current, b)
            CPlot.render()
        elif v == 2: # down
            b = CPlot.addRender2Zone(out[current], color='White')
            CPlot.replace(t, 1, current, b)
            current += 1
            if current >= nitems: current = 0
            b = CPlot.addRender2Zone(out[current], color='Red')
            CPlot.replace(t, 1, current, b)
            CPlot.render()
        elif v == 3:
            b = CPlot.addRender2Zone(out[current], color='White')
            CPlot.replace(t, 1, current, b)
            current -= 1
            if current < 0: current = nitems-1
            b = CPlot.addRender2Zone(out[current], color='Red')
            CPlot.replace(t, 1, current, b)
            CPlot.render()
        elif v == 4:
            b = CPlot.addRender2Zone(out[current], color='White')
            CPlot.replace(t, 1, current, b)
            current += 1
            if current >= nitems: current = 0
            b = CPlot.addRender2Zone(out[current], color='Red')
            CPlot.replace(t, 1, current, b)
            CPlot.render()
        elif v == 13:
            print('returning %s.'%str(current))
            return current, items[current]
        time.sleep(dt)
        l = ''; v = -1; CPlot.resetKeyboard()
