# - octree (pyTree) -
import Converter.Internal as Internal
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D
import Post.PyTree as P
import Transform.PyTree as T
import Converter.Mpi as Cmpi
import Apps.Coda.ToolboxIBM_CODA as TBX

import Connector.PyTree as X
rank = Cmpi.rank
NP = Cmpi.size
s = D.circle((0,0,0), 1., N=100)
o = G.octree([s], [0.8], balancing=1,dfar=5.)
a = T.splitNParts(o,N=NP, recoverBC=False)
for z in Internal.getZones(a):
    z[0]='ZONE_%d'%rank
    Cmpi._setProc(z,rank)

t = C.newPyTree(['Base']); t[2][1][2]=[Internal.getZones(a)[rank]]
tbb = Cmpi.createBBoxTree(t)
interDict = X.getIntersectingDomains(tbb)
graph = Cmpi.computeGraph(tbb, type='bbox', intersectionsDict=interDict, reduction=False)
a = Internal.getZones(t)[0];#a[0]='ZONE_%d'%rank

indicesFacesOrig=[]
extFaces = P.exteriorFaces(a, indicesFacesOrig)
extFaces[0]='extFaces_%d'%rank
indicesFacesOrig=indicesFacesOrig[0]
# Local hanging nodes
res = TBX.getHangingNodesInfoPara(a,extFaces, indicesFacesOrig, extFaces, indicesFacesOrig)
dictOfHangingNodes={}
if res[0] != []:
    dictOfHangingNodes[rank]=res

if NP>1:
    # Send info to opposite procs
    datas={}
    for opprank in graph[rank]:
        if opprank not in datas:
            datas[opprank]=[[extFaces, indicesFacesOrig]]
        else:
            datas[opprank].append([extFaces,indicesFacesOrig])

    destDatas=Cmpi.sendRecv(datas,graph)
    for i in destDatas:
        for res in destDatas[i]:
            extFacesOpp=res[0]; indicesFacesOrigOpp=res[1]
            res = TBX.getHangingNodesInfoPara(a, extFaces, indicesFacesOrig, extFacesOpp, indicesFacesOrigOpp)
            if res[0] != []:
                dictOfHangingNodes[i]=res
print("Hanging node on coarse mpi rank:", rank, ": ", dictOfHangingNodes)
