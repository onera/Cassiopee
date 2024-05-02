# - isoLine (array) -
import Post as P
import Converter as C
import Generator as G
import KCore.test as test

def F(x, y): return x*x+y*y

# Test sur un array
a = G.cartTetra( (0,0,0), (1,1,1), (10,10,1))
a = C.initVars(a, 'field', F, ['x','y'])
isos = []
min = C.getMinValue(a, 'field')
max = C.getMaxValue(a, 'field')
for v in range(20):
    value = min + (max-min)/18.*v
    try:
        i = P.isoLine(a, 'field', value)
        if i != []: isos.append(i)
    except: pass
test.testA(isos, 1)

# Test sur une liste d'arrays
b = G.cartTetra((12,0,0), (1,1,1), (10,10,1))
b = C.initVars(b, 'field', F, ['x','y'])
isos = []
min = C.getMinValue([a,b], 'field')
max = C.getMaxValue([a,b], 'field')
for v in range(20):
    value = min + (max-min)/18.*v
    try:
        i = P.isoLine([a,b], 'field', value)
        if i != []: isos.append(i)
    except: pass
test.testA(isos, 2)

# Essai sur des BARs
a = G.cartTetra((0,0,0), (1,1,1), (10,1,1))
a = C.initVars(a, 'field', F, ['x','y'])
isos = []
min = C.getMinValue(a, 'field')
max = C.getMaxValue(a, 'field')
for v in range(20):
    value = min + (max-min)/18.*v
    try:
        i = P.isoLine(a, 'field', value)
        if i != []: isos.append(i)
    except: pass
test.testA(isos, 3)
