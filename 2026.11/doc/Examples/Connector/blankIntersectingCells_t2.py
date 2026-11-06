# - blankIntersectingCells (array)
import Converter as C
import Generator as G
import Geom as D
import Transform as T
import Connector as X
import KCore.test as test
#-------------------------------------------------------------
# extrusion vers l exterieur : pas d intersection des facettes
#-------------------------------------------------------------
d = C.array('d', 5, 1, 1)
for i in range(1,d[2]+1): d[1][0,i-1] = 0.01*i
# structure
s = D.sphere6((0,0,0), 1.,N=10)
A = []; Ac = []
for i in range(len(s)):
    a = G.addNormalLayers(s[i], d)
    ac = C.node2Center(a); ac = C.initVars(ac,'cellN',1)
    Ac.append(ac); A.append(a)
res = X.blankIntersectingCells(A, Ac, tol=1.e-8)
test.testA(res,1)

# non structure hexa
sh = C.convertArray2Hexa(s)
A = []; Ac = []
for i in range(len(sh)):
    a = G.addNormalLayers(sh[i], d)
    ac = C.node2Center(a); ac = C.initVars(ac,'cellN',1)
    Ac.append(ac); A.append(a)
res = X.blankIntersectingCells(A, Ac, tol=1.e-8)
test.testA(res,2)

# non structure penta
st = C.convertArray2Tetra(s)
A = []; Ac = []
for i in range(len(st)):
    a = G.addNormalLayers(st[i], d)
    ac = C.node2Center(a); ac = C.initVars(ac,'cellN',1)
    Ac.append(ac); A.append(a)
res = X.blankIntersectingCells(A, Ac, tol=1.e-8)
test.testA(res,3)

#------------------------------------------
# extrusion vers l interieur : intersection
#------------------------------------------
d = C.array('d', 4, 1, 1)
for i in range(1,d[2]+1):
    d[1][0,i-1] =-0.1*i

# structure
s = D.sphere6((0,0,0), 1.,10)
A = []; Ac = []
for i in range(len(s)):
    a = G.addNormalLayers(s[i], d)
    ac = C.node2Center(a); ac = C.initVars(ac,'cellN',1)
    Ac.append(ac); A.append(a)
res = X.blankIntersectingCells(A, Ac, tol=1.e-8)
test.testA(res,4)

# non structure hexa
sh = C.convertArray2Hexa(s)
A = []; Ac = []
for i in range(len(s)):
    a = G.addNormalLayers(s[i], d)
    ac = C.node2Center(a); ac = C.initVars(ac,'cellN',1)
    Ac.append(ac); A.append(a)
res = X.blankIntersectingCells(A, Ac, tol=1.e-8)
test.testA(res,5)

# non structure tetra
st = C.convertArray2Tetra(s)
A = []; Ac = []
for i in range(len(s)):
    a = G.addNormalLayers(s[i], d)
    ac = C.node2Center(a); ac = C.initVars(ac,'cellN',1)
    Ac.append(ac); A.append(a)
res = X.blankIntersectingCells(A, Ac, tol=1.e-8)
test.testA(res,6)
