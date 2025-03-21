# - computeFrictionVelocity (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Connector.IBM as X_IBM
import numpy
import KCore.test as test

nb_node = 3

zsize = numpy.empty((1,3), Internal.E_NpyInt, order='F')
zsize[0,0] = nb_node; zsize[0,1] = 0; zsize[0,2] = 0
z = Internal.newZone(name='IBW_Wall',zsize=zsize,ztype='Unstructured')
coord_node = Internal.newGridCoordinates(parent=z)
# Image Points : PI
coord_node[2].append(['CoordinateX',numpy.array([0.1,0.2,0.3]),[],'DataArray_t'])
coord_node[2].append(['CoordinateY',numpy.array([0.1,0.2,0.3]),[],'DataArray_t'])
coord_node[2].append(['CoordinateZ',numpy.array([0.0005,0.0005,0.0005]),[],'DataArray_t'])

FS = Internal.newFlowSolution(parent=z)
# Image Points : PI
FS[2].append(['CoordinateX_PI',numpy.array([0.1,0.2,0.3]),[],'DataArray_t'])
FS[2].append(['CoordinateY_PI',numpy.array([0.1,0.2,0.3]),[],'DataArray_t'])
FS[2].append(['CoordinateZ_PI',numpy.array([0.0005,0.0005,0.0005]),[],'DataArray_t'])
# Wall Points : PW
FS[2].append(['CoordinateX_PW',numpy.array([0.1,0.2,0.3]),[],'DataArray_t'])
FS[2].append(['CoordinateY_PW',numpy.array([0.,0.,0.]),[],'DataArray_t'])
FS[2].append(['CoordinateZ_PW',numpy.array([0.0005,0.0005,0.0005]),[],'DataArray_t'])
# Quantities at Image Points
FS[2].append(['Density',numpy.array([1.,1.,1.]),[],'DataArray_t'])
FS[2].append(['ViscosityMolecular',numpy.array([5.e-8,5.e-8,5.e-8]),[],'DataArray_t'])
FS[2].append(['VelocityX',numpy.array([0.15,0.25,0.35]),[],'DataArray_t'])
FS[2].append(['VelocityY',numpy.array([0.,0.,0.]),[],'DataArray_t'])
FS[2].append(['VelocityZ',numpy.array([0.,0.,0.]),[],'DataArray_t'])
FS[2].append(['Pressure',numpy.array([1.,1.,1.]),[],'DataArray_t'])
FS[2].append(['utau',numpy.array([0.,0.,0.]),[], 'DataArray_t'])
FS[2].append(['yplus',numpy.array([0.,0.,0.]),[], 'DataArray_t'])

X_IBM._computeFrictionVelocity(z)
test.testT(z,1)
