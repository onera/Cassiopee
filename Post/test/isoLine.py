# - isoLine (array) -
import Post as P
import Converter as C
import Generator as G

def F(x, y): return x*x+y*y

a = G.cartTetra((0,0,0), (1,1,1), (10,10,1))
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
C.convertArrays2File([a]+isos, 'out.plt')
