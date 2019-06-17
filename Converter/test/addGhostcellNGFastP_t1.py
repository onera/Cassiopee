# - compute (pyTree) -
# Lamb vortex
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.PyTree as C
import Converter.Internal as Internal
import Connector.PyTree as X
import Converter.GhostCells as GC
import KCore.test as test

import FastP.PyTree as FastP
import Fast.PyTree as Fast
import Initiator.PyTree as I
import Connector.PyTree as X
import Post.PyTree as P
import CPlot.PyTree as CPlot
import sys

a = G.cartNGon((0,0,0), (1,1,1), (100,100,3))
t = C.newPyTree(['Base',a])
C._fillEmptyBCWith(t, 'extrap', 'BCExtrapolate', dim=3)

t = T.splitNParts(t, 2)

t = X.connectMatch(t,dim=3)


Internal._adaptNFace2PE(t, remove=False) 

#
ncouche =2
t = GC.addGhostCellsNG(t, nlayers=ncouche) 


test.testT(t, 1)

