# prepareAMRData : IBM preprocessing for HEXA/AMR grids
import Converter.PyTree as C
import Generator.PyTree as G
import Geom.PyTree as D
import Geom.IBM as D_IBM
import Connector.PyTree as X
import Transform.PyTree as T
a = D.naca(12.)
dimPb = 2
D_IBM._setSnear(a, 0.1)
D_IBM._setIBCType(a,"Musker")
D_IBM._setDfar(a, 20.)

tb = C.newPyTree(["BODY",a])
o = G.octree(tb, dfar=20, snearList=[0.1], balancing=1)
T._addkplane(o)
t = C.newPyTree(["CARTESIAN",o])
## ===========================================
##                 Prep IBM
## ===========================================
Reynolds     = 1.e6
yplus_target = 300
Lref         = 1.0
depth_IP     = 1
distance_IP  = D_IBM.computeModelisationHeight(Re=Reynolds, yplus=yplus_target, L=Lref)
print("distance_IP=", distance_IP,flush=True)
#
IBM_parameters = {
    "donor points":
    {
        "front type": "1",
        "depth DonorPoints": 1,  # for type 1
        # "distance DonorPoints": distance_DP,  # for type 2
    },
    "integration points":
    {
        "front type": "2",
        #"depth IntegrationPoints": depth_IP,  # for type 1
        "distance IntegrationPoints": distance_IP,  # for type 2
    },
    "IBM type":
    {
        "type": "global",
        "size elements body": 0.0050048828125,
    },
    "spatial discretization":
    {
        "type": "FV",
    }
}
t = X.prepareAMRData(tb, t, IBM_parameters=IBM_parameters, dim=dimPb)
C.convertPyTree2File(t,"out.cgns")
