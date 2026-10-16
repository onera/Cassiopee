# - generateAMRMesh (pyTree) -
import Converter.PyTree as C
import Generator.AMR as G_AMR
import Geom.PyTree as D
import Geom.IBM as D_IBM
import KCore.test as test

LOCAL = test.getLocal()

tb = C.convertFile2PyTree('tSymExtrude.cgns')
dictExtension = {
    "midPlane": -1
}
dictTolerance= {
    "midPlane": 1e-05
}
G_AMR._addExtensionInfo(tb, dictExtension, dictTolerance=dictTolerance)
test.testT(tb,1)

tbExtSym = G_AMR.createExtension__(tb)
test.testT(tbExtSym,2)

#C.convertPyTree2File(tbExtSym, 'check_tbExtSym.cgns')
