# - distance2Walls (array) -
# test : array structure, bodies structures, celln avec des 0 et des 2
import Dist2Walls
import Generator as G
import KCore.test as test
import Geom as D
import Converter as C

a = G.cart((0.,0.,0.),(0.1,0.1,0.1),(11,11,11))
sphere = D.sphere((0.5,0.5,0.5),0.2,30)
# pas de points masques
cellns = C.initVars(sphere, 'cellnf', 1.)
c = 1
for signed in [0,1]:
    for loc in ['nodes','centers']:
        for type in ['mininterf','ortho']:
            dist = Dist2Walls.distance2Walls(a, [sphere], cellnbodies=[cellns], loc=loc,type=type,signed=signed, dim=3)
            test.testA([dist],c)
            c += 1
# la moitie de la sphere est masquee
def cellN0__(x):
    if x < 0.5: return 0.
    else: return 1.
cellns = C.initVars(cellns,'cellN',cellN0__,['x'])

for signed in [0,1]:
    for loc in ['nodes','centers']:
        for type in ['mininterf','ortho']:
            dist = Dist2Walls.distance2Walls(a, [sphere], cellnbodies=[cellns], loc=loc,type=type,signed=signed, dim=3)
            test.testA([dist],c)
            c+=1

# la moitie de la sphere est interpolee
def cellN2__(x):
    if x < 0.5: return 2.
    else: return 1.
cellns = C.initVars(cellns,'cellN',cellN2__,['x'])

for signed in [0,1]:
    for loc in ['nodes','centers']:
        for type in ['mininterf','ortho']:
            dist = Dist2Walls.distance2Walls(a, [sphere], cellnbodies=[cellns], loc=loc, type=type, signed=signed, dim=3)
            test.testA([dist],c)
            c+=1
