# - Fast.IBM -
import Apps.Mesh.Cart as App
import Converter.Mpi as Cmpi

tloc = App.generateCartMesh('naca1DEuler.cgns', snears=0.01, dfar=10., dfarList=[], vmin=15, check=False, tbox=None, snearsf=None, ext=0, dimPb=2,sizeMax=1000)
Cmpi.convertPyTree2File(tloc,"cart.cgns")
