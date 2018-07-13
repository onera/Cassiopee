# test maskXRay : points de percage (pyTree)
import Connector.PyTree as X
import Geom.PyTree as D
import Transform.PyTree as T
import Converter.PyTree as C

# retourne les pierce pts sous forme de 'NODE'
surf = D.sphere((0,0,0), 0.5, 20)
surf = T.rotate(surf,(0.,0.,0.),(0.,1.,0.),90.)
xray = X.maskXRay__(surf)
