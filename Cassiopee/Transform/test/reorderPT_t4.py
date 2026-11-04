# - reorder (pyTree) -
# reorder une zone non carree avec topTree : les BCMatch doivent etre
# identiques a celles de connectMatch sur l'arbre reordonne
import Connector.PyTree as X
import Converter.Internal as Internal
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test
import Transform.PyTree as T


def records(t):
    out = {}
    for z in Internal.getZones(t):
        for g in Internal.getNodesFromType(z, "GridConnectivity1to1_t"):
            out[(z[0], Internal.getValue(g))] = [
                Internal.getNodeFromName1(g, n)[1].ravel().tolist()
                for n in ("PointRange", "PointRangeDonor", "Transform")
            ]
    return out


def case(dim, n, order, no):
    a = G.cart((0, 0, 0), (1, 1, 1), n)
    a[0] = "a"
    b = G.cart((n[0] - 1, 0, 0), (1, 1, 1), n)
    b[0] = "b"
    t = C.newPyTree(["Base", dim, a, b])
    t = X.connectMatch(t, dim=dim)
    T._reorder(Internal.getNodeFromName2(t, "a"), order, topTree=t)
    r = Internal.copyTree(t)
    Internal._rmNodesByType(r, "ZoneGridConnectivity_t")
    r = X.connectMatch(r, dim=dim)
    if records(t) != records(r):
        print(f"FAILED: reorder {order} on {n}")
    test.testT(t, no)


# 2D, rotation de 90 degres
case(2, (5, 9, 1), (2, -1, 3), 1)
# 3D, rotations
case(3, (5, 7, 9), (3, -1, -2), 2)
case(3, (5, 7, 9), (1, -3, 2), 3)
