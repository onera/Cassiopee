# - polyC1Mesher (array) -
import Converter as C
import Generator.PolyC1 as GP
import Generator as G
import Geom as D
import Transform as T
import KCore.test as test

# Definition de la geometrie
a = D.naca(12., 101)
a = T.reorder(a, (-1,2,3))

h = 0.1; hp = 0.001; density = 100.
m = GP.polyC1Mesher(a, h, hp, density)

for i in m[0]:
    v = G.getVolumeMap(i)
    min = C.getMinValue(v, 'vol')
    if (min <= 0):
        print('negative volume detected.')
test.testA(m[0], 1)

# avec raccord coincident : critere de split de la courbure faible
m = GP.polyC1Mesher(a, h, hp, density,splitCrit=0.01)
for i in m[0]:
    v = G.getVolumeMap(i)
    min = C.getMinValue(v, 'vol')
    if (min <= 0):
        print('negative volume detected.')
test.testA(m[0], 2)
