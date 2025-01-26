"""Interface with gmsh python module"""
try: import gmsh
except ImportError: raise ImportError("Gmsh: need gmsh module.")
import numpy

# Extrait les elements du maillage de dimension donnee et correspondant au tag donne
# si dim=-1 et tag=-1 extrait tout
def convertMesh2Arrays(dim=2, tag=-1):
    """Convert a gmsh mesh into an array."""
    model = gmsh.model
    # Recupere tous les noeuds du modele (toutes dimensions et tous les tags)
    # renvoie le tag des noeuds, les coords des noeuds et les param des noeuds
    (nodeTags, coords, params) = model.mesh.getNodes(dim=-1, tag=-1)
    # Recupere les connectivites des maillages de surface dans le modele
    # renvoie une liste de type d'elements, une liste de tag des elements et une liste
    # des tag des noeuds pour chaque element (connectivite)
    (elementTypes, elementTags, nodeElementTags) = model.mesh.getElements(dim=dim, tag=tag)

    # transforme coords en numpy pour Cassiopee
    npts = coords.size//3
    print('npts=',npts, coords.dtype)
    crds = numpy.empty((3,npts), dtype=numpy.float64)
    crds[0,0:npts] = coords[0:3*npts:3]
    crds[1,0:npts] = coords[1:3*npts:3]
    crds[2,0:npts] = coords[2:3*npts:3]

    # inverse nodeTags (bottleneck)
    inverse = {}
    for c, tag in enumerate(nodeTags): inverse[tag] = c+1

    # transforme elements en connectivite par elements pour Cassiopee
    out = []
    for c, ctype in enumerate(elementTypes):

        print('elt type=', ctype)
        prop = model.mesh.getElementProperties(elementType=ctype)
        eltName = prop[0]
        eltDim = prop[1]
        eltOrder = prop[2]
        nof = prop[3]

        if eltName == 'Line 1': eltType = 'BAR'
        elif eltName == 'Triangle 3': eltType = 'TRI'
        elif eltName == 'Tetrahedron 4': eltType = 'TETRA'
        elif eltName == 'Quadrilateral 4': eltType = 'QUAD'
        else:
            raise ValueError('convertMesh2Arrays: unknown element %s'%eltName)

        eltTags = elementTags[c]
        nodeEltTags = nodeElementTags[c]
        # nbre d'elements
        ne = eltTags.size
        kelts = numpy.empty( (nof,ne), dtype=numpy.int32)
        for n in range(nof): # bottleneck
            for i in range(ne):
                kelts[n,i] = inverse.get(nodeEltTags[nof*i+n], 1)

        a = ['x,y,z', crds, kelts, eltType]
        out.append(a)
    return out