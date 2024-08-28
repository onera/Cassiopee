# - volumeFromCrossSections (array) -
import Converter as C
import Geom as D
contours = []
for z in [0.,1.]:
    contours.append(D.circle((0,0,z),1.,N=15))
vol = D.volumeFromCrossSections(contours)
C.convertArrays2File([vol]+contours, 'out.plt')
