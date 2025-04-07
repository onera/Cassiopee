# - prepareIBMData Serial (pyTree) -
import Converter.PyTree as C
import Connector.IBM as X_IBM
import KCore.test as test

tb = C.convertFile2PyTree('naca1DNS.cgns')

# Prepare
vmin      = 42
dfars     = 5
snears    = 1
t, tc = X_IBM.prepareIBMData(tb             , None         , None     ,
                             snears=snears  , dfars=dfars  , vmin=vmin,
                             check=False    , frontType=1  , cartesian=False)

####
# The following lines are to avoid regression since the bug fix for duplicate information in tc
####
import Converter.Internal as Internal
for b in Internal.getBases(tc):
    for z in Internal.getZones(b):
        pos = 0
        z2 = Internal.copyRef(z)
        for zs in z2[2]:
            if 'ID' in zs[0] or 'IBCD' in zs[0]:
                Internal.addChild(z, zs, pos)
                pos +=2
            else:
                pos += 1
####

test.testT(t , 1)
test.testT(tc, 2)
#C.convertPyTree2File(t,'t_check.cgns')
