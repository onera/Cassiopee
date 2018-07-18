# - snapSharpEdges (array) -
import Generator as G
import Geom as D
import Transform as T
import KCore.test as test

# polylignes avec angles vifs
s = D.polyline([(0.02,0,0),(1,1,0),(2,1,0),(0.02,0,0)])
s = T.addkplane(s)

test.stdTestA(G.snapSharpEdges, [s], 1.)
