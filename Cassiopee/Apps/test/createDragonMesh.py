# Dragon mesh
import Converter.PyTree as C
import Apps.Mesh.Dragon as AppMesh
import Geom.PyTree as D

s = D.sphere((0,0,0),1,N=30)
dictOfParams={}
dictOfParams['hWall']=1e-3 # hauteur premiere maille
dictOfParams['h_factor']=4
dictOfParams['h_ext_factor']=4 #facteur*hauteur de la couche externe de prismes, -1

t = AppMesh.createDragonMesh(s, dictOfParams)
C.convertPyTree2File(t, "out.cgns")


#dictOfParams['sym_plane_axis']='Z'
#dictOfParams['sym_plane_values']=[0.]
#dictOfParams['niter']=100 # nb d iterations de lissage pour la couche de prismes
#dictOfParams['q']=raison # raison pour l extrusion des prismes, ex 1.2
#dictOfParams['nlayers']=nblayersforced # on force le nb de couches de prismes
#dictOfParams['h_ext']=snear # hauteur de la couche externe de prismes, -1 auto
#dictOfParams['dfar']=dfar # extension du domaine de calcul, auto si -1
#dictOfParams['NParts']=0 # decoupage en N parties  de la zone
#dictOfParams['holeFactor']=1.25 # facteur sur la hauteur de l evidement
