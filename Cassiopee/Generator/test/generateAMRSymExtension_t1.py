# - generateAMRMesh (pyTree) -
import Converter.PyTree as C
import Generator.AMR as G_AMR
import Geom.PyTree as D
import Geom.IBM as D_IBM
import KCore.test as test

LOCAL = test.getLocal()

tb = C.convertFile2PyTree('tSymExtrude.cgns')
tbv2 = G_AMR.createExtension__(tb)
test.testT(t,1)
