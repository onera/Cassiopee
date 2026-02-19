# - createDragonMesh (pyTree) -
import Apps.Mesh.Dragon as AppMesh
import Geom.PyTree as D
import KCore.test as test

s = D.sphere((0,0,0), 1, N=10)
dictOfParams={}
dictOfParams['hWall']=1.e-2 # hauteur premiere maille
dictOfParams['h_factor']=4
dictOfParams["niter"]=0
dictOfParams['h_ext_factor']=4 #fact*hauteur de la couche externe de prismes, -1
t = AppMesh.createDragonMesh(s, dictOfParams)
test.testT(t, 1)
