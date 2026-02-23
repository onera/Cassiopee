# - computeDiff (array) -
# ne marche pas car les connectivite + champ en noeud ne fonctionne pas
import Converter as C
import Post as P
import Generator as G
import KCore.test as test

def F(x):
    if x > 5.: return True
    else: return False
#
def celln(y):
    if y > 5.: return True
    else: return False

# computeDiff_t3 does not work for unstructured arrays, located at centers.
# #--------------
# # TRI
# #--------------
# ni = 30; nj = 40; nk = 1
# m = G.cartTetra((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
# m = C.node2Center(m)
# m = C.initVars(m, 'ro', F, ['x'])
# p = P.computeDiff(m,'ro')
# test.testA([p],1)
# #
# m = C.initVars(m, 'cellN', celln, ['y'])
# p = P.computeDiff(m,'ro')
# test.testA([p],2)

# #--------------
# # TETRA
# #--------------
# ni = 30; nj = 40; nk = 11
# m = G.cartTetra((0,0,0), (10./(ni-1),10./(nj-1),10./(nk-1)), (ni,nj,nk))
# m = C.initVars(m, 'ro', F, ['x'])
# p = P.computeDiff(m,'ro')
# test.testA([p],3)
# #
# m = C.initVars(m, 'cellN', celln, ['y'])
# p = P.computeDiff(m,'ro')
# test.testA([p],4)
# #--------------
# # QUAD
# #--------------
# ni = 30; nj = 40; nk = 1
# m = G.cartHexa((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
# m = C.initVars(m, 'ro', F, ['x'])
# p = P.computeDiff(m,'ro')
# test.testA([p],5)
# #
# m = C.initVars(m, 'cellN', celln, ['y'])
# p = P.computeDiff(m,'ro')
# test.testA([p],6)
# #--------------
# # HEXA
# #--------------
# ni = 30; nj = 40; nk = 11
# m = G.cartHexa((0,0,0), (10./(ni-1),10./(nj-1),10./(nk-1)), (ni,nj,nk))
# m = C.initVars(m, 'ro', F, ['x'])
# p = P.computeDiff(m,'ro')
# test.testA([p],7)
# #
# m = C.initVars(m, 'cellN', celln, ['y'])
# p = P.computeDiff(m,'ro')
# test.testA([p],8)
# #
# # PENTA
# #
# ni = 30; nj = 40; nk = 11
# m = G.cartPenta((0,0,0), (10./(ni-1),10./(nj-1),10./(nk-1)), (ni,nj,nk))
# m = C.initVars(m, 'ro', F, ['x'])
# p = P.computeDiff(m,'ro')
# test.testA([p],9)
# #
# m = C.initVars(m, 'cellN', celln, ['y'])
# p = P.computeDiff(m,'ro')
# test.testA([p],10)
# #
