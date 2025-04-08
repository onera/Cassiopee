# - getCurvatureRadius (array) -
import Geom as D
import Converter as C
pts = D.polyline([(6,0.01,1), (5.4,0.036,1), (4.8,0.064,1), (2.5,0.21,1),
                  (0.3,0.26,1),(0,0.047,1),(0,0,0)])
a = D.bezier(pts, 100)
rad = D.getCurvatureRadius(a)
C.convertArrays2File(a, 'out.plt')
