# - generateAMRMeshCartIn (pyTree) -
import Converter.PyTree as C
import Generator.AMR as G_AMR
import Geom.PyTree as D
import Geom.IBM as D_IBM
import KCore.test as test

LOCAL = test.getLocal()

vmins = [8,7,6,5]

# 2D
a = D.naca(12.)
dim = 2
D_IBM._setSnear(a, 1.e-03)
D_IBM._setIBCType(a,"Musker")
tb = C.newPyTree(["BODY",a])

dictGridCart={
    'gridType': 'cartesian',
    'cartbgExtent': [-20, -10, -20, 20, 10, 20],
    'cartbgBC': ['BCFarfield', 'BCFarfield', 'BCFarfield', 'BCFarfield', 'BCFarfield', 'BCFarfield'],
    'matchExtent': [False, False, False, True, True, True]
}

o, levelMax = G_AMR.generateCartBackgroundGrid(tb=tb, dim=dim, dictGridCart=dictGridCart)
test.testT(o, 1)
t = G_AMR.generateAMRMesh(tb=tb, levelMax=levelMax, vmins=vmins, dim=dim, tIn=o, localDir=LOCAL)
test.testT(t, 2)
#C.convertPyTree2File(t,'check_t.cgns')
