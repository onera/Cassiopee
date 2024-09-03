# - prepareIBMData MPI (pyTree) -
import Converter.PyTree as C
import Converter.Mpi as Cmpi
import Connector.IBM as X_IBM
import KCore.test as test

tb = C.convertFile2PyTree('../../Apps/test/naca1DNS.cgns')

# Prepare
vmin      = 42
dfars     = 5
snears    = 1
t, tc = X_IBM.prepareIBMData(tb             , None         , None     ,
                             snears=snears  , dfars=dfars  , vmin=vmin, 
                             check=False    , frontType=1  , cartesian=False)
if Cmpi.rank==0:
    test.testT(t , 1)
    test.testT(tc, 2)
