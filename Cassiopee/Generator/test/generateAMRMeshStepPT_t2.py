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

dimPb = 3
tb = C.convertFile2PyTree('m6cleanDouble.cgns')
tbox = C.convertFile2PyTree('m6tbox.cgns')

##
isRunMeshGen = False
dim = dimPb
check = False
levelMax = 0
NumMinDxLarge = 2
octreeMode = 1
vmins = [10,8,6]
snears = 2e-02
vminsTbox = 3
dfars = 20


D_IBM._setDfar(tb,dfars)
for bNum, b in enumerate(Internal.getBases(tb)):
    D_IBM._setSnear(b, snears*(bNum+1))
D_IBM._symmetrizePb(tb, ['Base', 'Base.0'], snears, 2)
D_IBM._setSnear(tbox, 8*snears)

#Note: the 3 next checks are inside G_AMR.checkInputIBM__
tb_noSym, tbExtSym, dir_sym = G_AMR.checkBaseSym__(tb, None, dim, check, LOCAL)
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
o, tCartIn, levelMax, snears, dfars, vmins = G_AMR.generateCoarseGrid__(tb, tb_tbox_noSym, None, snears, snearsTb, nbases, vmins, dfars,
                                                                        dim, dir_sym, levelMax, octreeMode, NumMinDxLarge)
if Cmpi.master:
    test.testT(o, 10)
    test.testO(vmins, 11)
    test.testO(snears, 12)
    test.testO(dfars, 13)

if isRunMeshGen:
    if Cmpi.master: C.convertPyTree2File(o, LOCAL+'/skeleton.cgns')
    Cmpi.barrier()
   
if isRunMeshGen:
    toffset = G_AMR.generateOffset__(tb_tbox, nbases, dfars, vmins, snears, dim, False, nboxes, tbExtSym, 'xray', tCartIn)
    Cmpi.barrier()
    hmin_skel = (C.getMinValue(o,"centers:vol"))**(1/dim)
    hmin = hmin_skel * 2 ** (-levelMax)
    o = G_AMR.adaptMesh__(LOCAL+'/skeleton.cgns', hmin, tb_noSym, toffset=toffset, dim=dim, nboxes=nboxes)
    Cmpi.convertPyTree2File(o, 'check_t.cgns')

