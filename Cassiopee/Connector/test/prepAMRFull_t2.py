import Converter.PyTree as C
import Geom.IBM as D_IBM
import Connector.AMR as X_AMR
import KCore.test as test
LOCAL = test.getLocal()
LOCAL +='/'
## ===================
##       INPUTS
## ===================
snear      = 0.002
dFar       = 10
levelMax   = 20
vmins      = [25,24,24,20,21,5,4]
dim        = 2

################### IBM pre-processing #########################
IBM_parameters = {
    "donor points": {
        "front type": "1",
        "depth DonorPoints": 1,  # for type 1
    },
    "integration points": {
        "front type": "1",
        "depth IntegrationPoints": 1,  # for type 1
    },
    "IBM type": {
        "type": "global",
        "size elements body": snear,
    },
    "spatial discretization": {
        "type": "FV",
    },
}

tb = D_IBM.naca0012(snear, 'Musker', 0.)
D_IBM._setDfar(tb, dFar)
##C.convertPyTree2File(tb, 'tb.cgns')

tIBM = X_AMR.prepareAMRIBM(tb, levelMax, vmins, dim, IBM_parameters, OutputAMRMesh=False, check=False, localDir=LOCAL)
test.testT(tIBM,1)
