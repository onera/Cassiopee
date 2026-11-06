# - meshSize (pyTree) -
import Initiator.PyTree as I

hp = I.meshSize(UInf=1., RoInf=1.223, ReInf=6000, LInf=1., yplus=1, algo='Turbulent')
print(hp)
