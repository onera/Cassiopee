# - connect1D (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C

P1 = [-0.5,0,0]; P1b = [0.5,0,0]
P2 = [1,-1.5,0]; P2b = [1,-0.5,0]
l1 = D.line(P1,P1b)
l2 = D.line(P2,P2b)

out = D.connect1D([l1,l2], sharpness=0)
C.convertPyTree2File(out, 'out.cgns')
