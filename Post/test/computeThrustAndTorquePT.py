# - computeThrustAndTorque (pyTree) -
import Converter.PyTree as C
import Post.Rotor as PR
import math

Mtip = 0.6462; MU = 0.4
CHORD = 0.14; AR = 2.1
RoInf = 1.225; PInf = 101325.; AINF = 340.1
SIGMA = 4*CHORD / (math.pi*AR)
teff = C.convertFile2PyTree('stress_419.cgns')

thrust,torque = PR.computeThrustAndTorque(teff, 419., PInf, center=(0,0,0), relativeShaft=0.); print('thrust', thrust, 'torque',torque)
#> thrust [368.952736931205, -39.172151751517326, 3543.2002154791667] torque [226.5518638611441, 950.9479780913017, -935.2791149345967]

accu = {}
thrust, torque = PR.computeThrustAndTorque(teff, 419., PInf, center=(0,0,0), relativeShaft=0., accumulatorThrust=accu)
ret = PR.exportAccumulatorPerRadius(accu, vars=['ThrustX','ThrustY','ThrustZ','TorqueX','TorqueY','TorqueZ'])
C.convertPyTree2File(ret, 'Thrust.cgns')
