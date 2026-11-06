# - polyC1Mesher (array) -
import Converter as C
import Generator.PolyC1 as GP
import Generator as G
import Transform as T

# Read geometry from svg file
a = C.convertFile2Arrays('Data/curve.svg', density=1)[0]
a = T.homothety(a,(0,0,0),0.01)
a = T.reorder(a, (1,2,3))

h = 0.2; hp = 0.001; density = 10.; splitCrit = 2.
m = GP.polyC1Mesher(a, h, hp, density, splitCrit)

for i in m[0]:
    v = G.getVolumeMap(i)
    min = C.getMinValue(v, 'vol')
    if min <= 0:
        print('negative volume detected.')

C.convertArrays2File(m[0], 'out.plt')
