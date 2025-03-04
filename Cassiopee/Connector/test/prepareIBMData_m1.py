# - prepareIBMData MPI (pyTree) -
import Converter.PyTree as C
import Converter.Mpi as Cmpi
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
if Cmpi.rank==0:
    import Converter.Internal as Internal
    ####
    # The following lines are to avoid regression since the bug fix for duplicate information in tc
    ####
    dictOfDoublons = {
        'Cart.2X1':['ID_Cart.4X0', 'ID_CartX0', 'ID_Cart.3X1'],
        'Cart.3X1':['ID_Cart.4X0', 'ID_Cart.2X1'],
        'Cart.3X0':['ID_Cart.5X0', 'ID_Cart.1X0'],
        'Cart.6X0':['ID_Cart.0X0', 'ID_Cart.0X1', 'ID_Cart.5X0'],
        'Cart.1X1':['IBCD_3_Cart.1X1']
    }

    for b in Internal.getBases(tc):
        for z in Internal.getZones(b):
            if z[0] in dictOfDoublons:
                pos = 0
                z2 = Internal.copyRef(z)
                for zs in z2[2]:
                    if ('ID' in zs[0] or 'IBCD' in zs[0]) and zs[0] in dictOfDoublons[z[0]]:
                        Internal.addChild(z, zs, pos)
                        pos +=2
                    else:
                        pos += 1
    ####

    test.testT(t , 1)
    test.testT(tc, 2)
