# - prepareIBMData Serial (pyTree) -
import Converter.PyTree as C
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
test.testT(t , 1)
test.testT(tc, 2)
#C.convertPyTree2File(t,'t_check.cgns')
