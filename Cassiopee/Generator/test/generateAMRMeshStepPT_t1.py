# - generateAMRMesh (pyTree) -
import Converter.PyTree as C
import Converter.Mpi as Cmpi
import Converter.Internal as Internal
import Transform.PyTree as T
import Generator.AMR as G_AMR
import Generator.PyTree as G
import Geom.PyTree as D
import Geom.IBM as D_IBM
import KCore.test as test

LOCAL = test.getLocal()

# 2D
a = D.naca(12.)
dimPb = 2
D_IBM._setSnear(a, 0.001)
D_IBM._setIBCType(a,"Musker")
D_IBM._setDfar(a, 10.)

b = Internal.copyTree(a)
T._translate(b, (0.,-1.,0.))
D_IBM._setSnear(b, 0.004)

tb = C.newPyTree(["Base1", a, "Base2", b])

#Ref 1 & 2
a1 = D.line((1.5,-0.2,0), (2.5,-0.2,0))
a2 = D.line((2.5,-0.2,0), (2.5, 0.8,0))
a3 = D.line((2.5, 0.8,0), (1.5, 0.8,0))
a4 = D.line((1.5, 0.8,0), (1.5,-0.2,0))

a = T.join([a1,a2,a3,a4])
D_IBM._setDfar(a, 10.)
D_IBM._setSnear(a, 0.004)

b = Internal.copyTree(a)
T._translate(b, (0.,-1.5,0.))
D_IBM._setSnear(b, 0.008)

tbox = C.newPyTree(["tbox1", a, "tbox2", b])

##
isRunMeshGen = False
dim = dimPb
check = False
levelMax = 0
NumMinDxLarge = 2
octreeMode = 1
vmins = [[25,24,24,20,21,5,4]]
snears = 1e-03
vminsTbox = 3
dfars = 10

tb_noSym, tbExtSym, dir_sym = G_AMR.checkBaseSym(tb, None, dim, check, LOCAL)
if Cmpi.master:
    test.testT(tb_noSym,1)
snears, snearsTb, vmins, nbases  = G_AMR.checkSnearsVmins__(tb_noSym, tbox, snears, vmins, levelMax)
if Cmpi.master:
    test.testO(snears, 2)
    test.testO(snearsTb, 3)
    test.testO(vmins, 4)
tb_tbox, tb_tbox_noSym, vmins, snears, nbases, nboxes = G_AMR.checkTbox__ (tb, tb_noSym, tbox, vmins, snears, nbases, vminsTbox, levelMax)
if Cmpi.master:
    test.testT(tb_tbox, 5)
    test.testT(tb_tbox_noSym, 6)
    test.testO(snears, 7)
    test.testO(snearsTb, 8)
    test.testO(vmins, 9)
o, tCartIn, extrude, levelMax, vmins = G_AMR.generateCoarseGrid__(tb, None, snears, snearsTb, nbases, vmins, dfars,
                                                                dim, levelMax, octreeMode)
if Cmpi.master:
    test.testT(o, 10)
    test.testO(vmins, 11)

if isRunMeshGen:
    if Cmpi.master: C.convertPyTree2File(o, LOCAL+'/skeleton.cgns')
    Cmpi.barrier()

G._getVolumeMap(o)
snears = G_AMR.updateSnears__(o, snears, dim, levelMax, nbases)
dfars = G_AMR.updateDfars__(o, tb_tbox_noSym, dim, dir_sym, tCartIn, extrude, NumMinDxLarge)
if Cmpi.master:
    test.testO(snears, 12)
    test.testO(dfars, 13)
toffset = G_AMR.generateOffset__(tb_tbox, nbases, dfars, vmins, snears, dim, False, nboxes, tbExtSym, 'xray', tCartIn)
if Cmpi.master:
    test.testT(toffset, 14)

if isRunMeshGen:
    Cmpi.barrier()
    hmin_skel = (C.getMinValue(o,"centers:vol"))**(1/dim)
    hmin = hmin_skel * 2 ** (-levelMax)
    o = G_AMR.adaptMesh__(LOCAL+'/skeleton.cgns', hmin, tb_noSym, toffset=toffset, dim=dim, nboxes=nboxes)
    C.convertPyTree2File(o, 'check_t.cgns')

