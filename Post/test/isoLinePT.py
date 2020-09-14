# - isoLine (pyTree) -
import Post.PyTree as P
import Converter.PyTree as C
import Generator.PyTree as G

def F(x, y): return x*x+y*y

a = G.cartTetra( (0,0,0), (1,1,1), (10,10,1))
a = C.initVars(a, 'field', F, ['CoordinateX','CoordinateY'])

isos = []
min = C.getMinValue(a, 'field')
max = C.getMaxValue(a, 'field')
for v in range(20):
    value = min + (max-min)/18.*v
    try:
        i = P.isoLine(a, 'field', value)
        isos.append(i)
    except: pass
t = C.newPyTree(['Base',3,'ISOS',1])
t[2][1][2].append(a)
t[2][2][2] += isos
C.convertPyTree2File(t, 'out.cgns')
