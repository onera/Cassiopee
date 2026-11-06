# - optimizeOverlap (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T
import Connector.PyTree as X
import KCore.test as test

Ni = 50; Nj = 50; Nk = 2
a = G.cart((0,0,0),(1./(Ni-1), 1./(Nj-1),1), (Ni,Nj,Nk))
b = G.cart((0,0,0),(2./(Ni-1), 2./(Nj-1),1), (Ni,Nj,Nk)); b[0] = 'cart2'
a = T.rotate(a, (0,0,0), (0,0,1), 10.); a = T.translate(a, (0.5,0.5,0))
a = C.addBC2Zone(a,'wall','BCWall','imin')
b = C.addBC2Zone(b,'wall','BCWall','imin')

t = C.newPyTree(['Base1','Base2'])
t[2][1][2].append(a); t[2][2][2].append(b)
t[2][1] = C.addState(t[2][1], 'EquationDimension',2)
t[2][2] = C.addState(t[2][2], 'EquationDimension',2)
t = C.fillEmptyBCWith(t,'overlap','BCOverlap', dim=2)
t = C.addVars(t,'Density'); t = C.addVars(t,'centers:Pressure')
t2 = X.optimizeOverlap(t)
test.testT(t2)
#
# priorite sur Base2
#
t2 = X.optimizeOverlap(t, priorities=['Base1',1,'Base2',0])
test.testT(t2,2)

# Like the first test, but using an intersection dictionnary:
interDict = X.getIntersectingDomains(t, method='hybrid')
t2 = X.optimizeOverlap(t, intersectionsDict=interDict)
test.testT(t2,3)

# Like the second test, but using an Intersection Dictionnary:
t2 = X.optimizeOverlap(t, priorities=['Base1',1,'Base2',0],
                       intersectionsDict=interDict)
test.testT(t2,4)
