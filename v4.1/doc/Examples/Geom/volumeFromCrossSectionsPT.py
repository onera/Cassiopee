# - volumeFromCrossSection (pyTree) -
import Converter.PyTree as C
import Geom.PyTree as D
contours = []
for z in [0.,1.]:
    contours.append(D.circle((0,0,z),1.,N=15))
vol = D.volumeFromCrossSections(contours)
C.convertPyTree2File(vol,'out.cgns')
