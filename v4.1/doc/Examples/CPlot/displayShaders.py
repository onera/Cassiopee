# - display (pyTree) -
# Affichage des shaders (mode solid)
import Geom.PyTree as D
import CPlot.PyTree as CPlot
import Transform.PyTree as T
import Converter.PyTree as C
import time

a = D.sphere((0,0,0),1.)
a = CPlot.addRender2Zone(a, material='Chrome', color='White')

b = T.translate(a, (0,2.5,0))
b = CPlot.addRender2Zone(b, material='Glass', color='Blue')

c = T.translate(a, (0,5.,0))
c = CPlot.addRender2Zone(c, material='Wood', color='Brown')

d = T.translate(a, (0,0,-2.5))
d = CPlot.addRender2Zone(d, material='Marble', color='Yellow')

e = T.translate(a, (0,2.5,-2.5))
e = CPlot.addRender2Zone(e, material='XRay', color='Blue')

f = T.translate(a, (0,5,-2.5))
f = CPlot.addRender2Zone(f, material='Metal', color='Blue')

i = T.translate(a, (0,0,-5))
i = CPlot.addRender2Zone(i, material='Granite', color='White')

g = T.translate(a, (0,2.5,-5.))
g = C.initVars(g, '{F}={CoordinateX}*{CoordinateX}+{CoordinateY}*{CoordinateY}')
g = C.node2Center(g, 'F')
g = CPlot.addRender2Zone(g, material='Solid', color='Iso:F')

h = T.translate(a, (0,2.5,-5.))
h = C.initVars(h, '{F}={CoordinateX}*{CoordinateX}+{CoordinateY}*{CoordinateY}')
h = C.node2Center(h, 'F')
h = C.rmVars(h, ['F'])
h = T.translate(h, (0,2.5,0))
h = CPlot.addRender2Zone(h, material='Solid', color='Iso:centers:F')

t = C.newPyTree(['Base',a,b,c,d,e,f,g,h,i])

# Tests des shaders
CPlot.display(t, mode=2) # mode=2 active les shaders
time.sleep(2)

# Ajout du overlay mesh
CPlot._addRender2Zone(t, meshOverlay=1)
CPlot.display(t, mode=2)
time.sleep(2)

# Ajout du blending
CPlot._addRender2Zone(t, meshOverlay=0)
for i in range(30):
    CPlot._addRender2Zone(t, blending=1.-i*1./30.)
    CPlot.display(t, mode=2)
    time.sleep(0.1)
