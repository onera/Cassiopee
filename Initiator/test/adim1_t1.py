# - adim1 -
import KCore.test as test
import Initiator.Adim as Adim
state = Adim.adim1(MInf=0.8, alphaZ=1., ReInf=1.e6)
test.testO(state, 1)
