#include <string.h>
#include "post.h"

#define MAXNEIS 12

using namespace std;
using namespace K_FLD;

static E_Float **mallocF(E_Int rows, E_Int cols)
{
    E_Float **array = (E_Float **)malloc(rows * sizeof(E_Float *));
    for (E_Int i = 0; i < rows; i++)
        array[i] = (E_Float *)malloc(cols * sizeof(E_Float));
    return array;
}

static void matrixMultiply
(
    E_Float **A, E_Float **B, E_Float **res,
    E_Int rowA, E_Int colA, E_Int colB
)
{
    for (E_Int i = 0; i < rowA; i++) {
        for (E_Int j = 0; j < colB; j++) {
            res[i][j] = 0;
            for (E_Int k = 0; k < colA; k++)
                res[i][j] += A[i][k] * B[k][j];
        }
    }
}

static void matVecMultiply
(
    E_Float **A, E_Float *V, E_Float *res,
    E_Int rowA, E_Int colA
)
{
    for (E_Int i = 0; i < rowA; i++) {
        res[i] = 0;
        for (E_Int j = 0; j < colA; j++) {
            res[i] += A[i][j] * V[j]; 
        }
    }
}

static void matrix2Inverse(E_Float **A, E_Float **invA)
{
    E_Float a11 = A[0][0]; E_Float a12 = A[0][1];
    E_Float a21 = A[1][0]; E_Float a22 = A[1][1];

    E_Float det = a11*a22 - a12*a21;
    E_Float invDet = 1 / det;

    invA[0][0] = a22 * invDet;
    invA[0][1] = -a12 * invDet;
    invA[1][0] = -a21 * invDet;
    invA[1][1] = a11 * invDet;
}

static void matrix3Inverse(E_Float **A, E_Float **invA)
{
    E_Float a11 = A[0][0]; E_Float a12 = A[0][1]; E_Float a13 = A[0][2];
    E_Float a21 = A[1][0]; E_Float a22 = A[1][1]; E_Float a23 = A[1][2];
    E_Float a31 = A[2][0]; E_Float a32 = A[2][1]; E_Float a33 = A[2][2];

    E_Float det = a11*(a33*a22 - a32*a23)
                - a21*(a33*a12 - a32*a13)
                + a31*(a23*a12 - a22*a13);
    
    E_Float invDet = 1 / det;

    invA[0][0] = (a33*a22 - a32*a23) * invDet;
    invA[0][1] = (a32*a13 - a33*a12) * invDet;
    invA[0][2] = (a23*a12 - a22*a13) * invDet;
    invA[1][0] = (a31*a23 - a33*a21) * invDet;
    invA[1][1] = (a33*a11 - a31*a13) * invDet;
    invA[1][2] = (a21*a13 - a23*a11) * invDet;
    invA[2][0] = (a32*a21 - a31*a22) * invDet;
    invA[2][1] = (a31*a12 - a32*a11) * invDet;
    invA[2][2] = (a22*a11 - a21*a12) * invDet;
}

static void matrixInverse(E_Int dim, E_Float **A, E_Float **invA)
{
    /* Assumes A is invertible */
    if (dim == 2) matrix2Inverse(A, invA);
    else matrix3Inverse(A, invA);
}

PyObject *K_POST::computeGradLSQ(PyObject *self, PyObject *args)
{
    PyObject *array;
    PyObject *arrayc;
    E_Int *DIM = new E_Int;

    if (!PyArg_ParseTuple(args, "OOi", &array, &arrayc, DIM)) return NULL;

    char *varString;
    char *eltType;
    FldArrayF *f;
    FldArrayI *cn;
    E_Int ni, nj, nk;
    E_Int posx = -1, posy = -1, posz = -1;

    E_Int res = K_ARRAY::getFromArray(array, varString, f, ni, nj, nk, cn, eltType, true);

    if (res != 1 && res != 2) {
        PyErr_SetString(PyExc_TypeError, "computeGradLSQ: invalid array.");
        return NULL;
    }

    if (res == 1 || strcmp(eltType, "NGON") != 0) {
        RELEASESHAREDB(res, array, f, cn);
        PyErr_SetString(PyExc_TypeError, "computeGradLSQ: only for NGons.");
        return NULL;
    }

    posx = K_ARRAY::isCoordinateXPresent(varString);
    posy = K_ARRAY::isCoordinateYPresent(varString);
    posz = K_ARRAY::isCoordinateZPresent(varString);

    if (posx == -1 || posy == -1 || posz == -1) {
        PyErr_SetString(PyExc_TypeError, "computeGradLSQ: coordinates not found in array.");
        RELEASESHAREDB(res, array, f, cn);
        return NULL;
    }

    posx++; posy++; posz++;

    char *varStringc; char *eltTypec;
    FldArrayF *fc; FldArrayI *cnc;
    E_Int nic, njc, nkc;
    res = K_ARRAY::getFromArray(arrayc, varStringc, fc, nic, njc, nkc, cnc, eltTypec, true);

    E_Int npts = f->getSize();
    printf("npts: %d\n", npts);

    E_Int nfld = fc->getNfld();
    vector<char *> vars;
    K_ARRAY::extractVars(varStringc, vars);

    vector<char *> varStrings;
    for (E_Int i = 0; i < nfld; i++) {
        char *local;
        computeGradVarsString(vars[i], local);
        varStrings.push_back(local);
    }

    E_Int size = 0;
    for (E_Int i = 0; i < nfld; i++)
        size += strlen(varStrings[i]) + 1;
    
    char *varStringOut = new char[size];
    char *pt = varStringOut;
    for (E_Int i = 0; i < nfld; i++) {
        char *v = varStrings[i];
        for (E_Int j = 0; j < strlen(v); j++) {
            *pt = v[j];
            pt++;
        }
        *pt = ',';
        pt++;
        delete [] varStrings[i];
    }
    pt--;
    *pt = '\0';
    for (size_t i = 0; i < vars.size(); i++)
        delete [] vars[i];

    FldArrayI cFE;
    K_CONNECT::connectNG2FE(*cn, cFE);
    E_Int *cFE1 = cFE.begin(1);
    E_Int *cFE2 = cFE.begin(2);
    E_Int* cnp = cn->begin();
    E_Int nfaces = cnp[0];

    FldArrayI cNFace;
	E_Int nelts;
	K_CONNECT::connectFE2NFace(cFE, cNFace, nelts);

    printf("nfaces: %d\n", nfaces);
    printf("nelts: %d\n", nelts);

    vector<vector<E_Int>> cell2Cells(nelts);
    E_Int own, nei;
	for (E_Int i = 0; i < nfaces; i++) {
		own = cFE1[i]-1;
		nei = cFE2[i]-1;

		if (nei < 0) {
			continue;
		}
		
		cell2Cells[own].push_back(nei);
		cell2Cells[nei].push_back(own);
	}

	FldArrayI posFaces;
	K_CONNECT::getPosFaces(*cn, posFaces);

	E_Int pos, stride;
	E_Int point;
	E_Float *pX = f->begin(posx);
	E_Float *pY = f->begin(posy);
	E_Float *pZ = f->begin(posz);
	FldArrayF face_coords;
	face_coords.malloc(nfaces, 3);
	face_coords.setAllValuesAtNull();
    cnp = cn->begin();

	for (E_Int i = 0; i < nfaces; i++) {
		cnp = cn->begin();
		pos = posFaces[i];
		cnp += pos;
		stride = cnp[0];
		for (E_Int j = 1; j <= stride; j++) {
			point = cnp[j]-1;	
			face_coords(i, 1) += pX[point]; 
			face_coords(i, 2) += pY[point];
			face_coords(i, 3) += pZ[point];
		}
		face_coords(i, 1) /= (E_Float) stride;
		face_coords(i, 2) /= (E_Float) stride;
		face_coords(i, 3) /= (E_Float) stride;
	}

	FldArrayF cellCenters;
	cellCenters.malloc(nelts, 3);
	cellCenters.setAllValuesAtNull();
	FldArrayI posCells;
	K_CONNECT::getPosElts(*cn, posCells);
	E_Int face;
	for (E_Int i = 0; i < nelts; i++) {
		cnp = cn->begin();
		pos = posCells[i];
		cnp += pos;
		stride = cnp[0];
		for (E_Int j = 1; j <= stride; j++) {
			face = cnp[j]-1;
			cellCenters(i, 1) += face_coords(face, 1);
			cellCenters(i, 2) += face_coords(face, 2);
			cellCenters(i, 3) += face_coords(face, 3);
		}
		cellCenters(i, 1) /= (E_Float) stride;
		cellCenters(i, 2) /= (E_Float) stride;
		cellCenters(i, 3) /= (E_Float) stride;
	}

    PyObject *tpl = K_ARRAY::buildArray(nfld*3, varStringOut, npts, nelts, -1, eltType, true,
        cn->getSize() * cn->getNfld());
    E_Int *cnnp = K_ARRAY::getConnectPtr(tpl);
    K_KCORE::memcpy__(cnnp, cn->begin(), cn->getSize()*cn->getNfld());
    E_Float *gnp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF gp(nelts, 3*nfld, gnp, true);
    gp.setAllValuesAtNull();

    // Here we iterate over cells
    // Could be done iterating over faces: faster but needs more memory
    E_Float Grad[3] = {0, 0, 0};
    E_Float **d = mallocF(MAXNEIS, *DIM);
    E_Float **td = mallocF(*DIM, MAXNEIS);
    E_Float **LSQ = mallocF(*DIM, MAXNEIS);
    E_Float **G = mallocF(*DIM, *DIM);
    E_Float **invG = mallocF(*DIM, *DIM);
    E_Float *deltaFld = (E_Float *)malloc(MAXNEIS * sizeof(E_Float));
    
    E_Float cx, cy, cz;
    E_Int nNeis;
    E_Float *fcp;

    for (E_Int fld = 0; fld < nfld; fld++) {
        fcp = fc->begin(fld+1);
        E_Float* gpx = gp.begin(3*fld+1);
        E_Float* gpy = gp.begin(3*fld+2);
        E_Float* gpz = gp.begin(3*fld+3);

        for (E_Int i = 0; i < nelts; i++) {
            nNeis = cell2Cells[i].size();

            if (nNeis > MAXNEIS) {
                fprintf(stderr, "Cell %d exceeded allocated number of neighbours\n", i);
                return 0;
            }

            cx = cellCenters(i, 1);
            cy = cellCenters(i, 2);
            cz = cellCenters(i, 3);

            // d
            for (E_Int j = 0; j < nNeis; j++) {
                nei = cell2Cells[i][j];
                d[j][0] = cellCenters(nei, 1) - cx;
                d[j][1] = cellCenters(nei, 2) - cy;
                d[j][2] = cellCenters(nei, 3) - cz;
                deltaFld[j] = fcp[nei] - fcp[i];
            }

            // td
            for (E_Int j = 0; j < *DIM; j++) {
                for (E_Int k = 0; k < nNeis; k++) {
                    td[j][k] = d[k][j];
                }
            }

            // G = td . d
            matrixMultiply(td, d, G, *DIM, nNeis, *DIM);
            matrixInverse(*DIM, G, invG);

            // LSQ = invG . td
            matrixMultiply(invG, td, LSQ, *DIM, *DIM, nNeis);

            // Gradi = LSQ . deltaFld
            matVecMultiply(LSQ, deltaFld, Grad, *DIM, nNeis);

            gpx[i] = Grad[0];
            gpy[i] = Grad[1];
            gpz[i] = Grad[2];
        }
    }

    RELEASESHAREDB(res, array, f, cn);
    RELEASESHAREDB(res, arrayc, fc, cnc);

    delete [] varStringOut;
    for (E_Int i = 0; i < MAXNEIS; i++) free(d[i]);
    for (E_Int i = 0; i < *DIM; i++) {
        free(td[i]);
        free(LSQ[i]);
        free(G[i]);
        free(invG[i]);
    }
    free(d);
    free(td);
    free(LSQ);
    free(G);
    free(invG);
    free(deltaFld);
    free(DIM);

    return tpl;
}