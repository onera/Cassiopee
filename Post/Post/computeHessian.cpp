#include <string.h>
#include "post.h"
#include <unordered_set>

using namespace K_FLD;
using namespace std;

static void computeHessianVarsString(char* varString, char*& varStringOut)
{
	vector<char*> vars;
	K_ARRAY::extractVars(varString, vars);
	E_Int c = -1;
	E_Int varsSize = vars.size();
	E_Int sizeVarStringOut = 0;
	
	for (E_Int v = 0; v < varsSize; v++) {
		E_Int vsize = strlen(vars[v]);
   		sizeVarStringOut += vsize+4; // HxxvarString,
	}

  	varStringOut = new char [9 * sizeVarStringOut];

  	for (E_Int v = 0; v < varsSize; v++) {
		char*& var0 = vars[v];
		if (strcmp(var0, "x") != 0 && strcmp(var0, "y") != 0 && strcmp(var0, "z") != 0) {
      			if (c == -1) {
        			strcpy(varStringOut, "Hxx");
				c = 1;
			} else {
				strcat(varStringOut, ",Hxx");
			}
			strcat(varStringOut, var0);
			strcat(varStringOut, ",Hxy");
			strcat(varStringOut, var0);
			strcat(varStringOut, ",Hxz");
			strcat(varStringOut, var0);
			strcat(varStringOut, ",Hyx");
			strcat(varStringOut, var0);
			strcat(varStringOut, ",Hyy");
			strcat(varStringOut, var0);
			strcat(varStringOut, ",Hyz");
			strcat(varStringOut, var0);
			strcat(varStringOut, ",Hzx");
			strcat(varStringOut, var0);
			strcat(varStringOut, ",Hzy");
			strcat(varStringOut, var0);
			strcat(varStringOut, ",Hzz");
			strcat(varStringOut, var0);
		}
	}
	
	for (E_Int v = 0; v < varsSize; v++) delete [] vars[v];
}

#define MAXNEIS 14

// Calcul de la hessienne d'un ensembe de champs en centres
// La hessienne est founie aux centres des cellules

PyObject *K_POST::computeHessian(PyObject *self, PyObject *args)
{
	PyObject *array;
	PyObject *arrayc;
	E_Int *DIM = new E_Int;
	
	if (!PyArg_ParseTuple(args, "OOi", &array, &arrayc, DIM)) {
		printf("Error in parsing args\n");
		return NULL;
	}

	if (*DIM != 2) {
		PyErr_SetString(PyExc_TypeError, "computeHessian: only for 2D.");
		return NULL;
	}

	// Check array
	char *varString;
	char *eltType;
	FldArrayF *f;
	FldArrayI *cn;
	E_Int ni, nj, nk;
	E_Int posx = -1;
	E_Int posy = -1;
	E_Int posz = -1;
	E_Int res = K_ARRAY::getFromArray(array, varString, f, ni, nj, nk, cn, eltType, true);

	if (res != 1 && res != 2) {
		PyErr_SetString(PyExc_TypeError, "computeHessian: invalid array.");
		return NULL;
	}

	if (res == 1 || strcmp(eltType, "NGON") != 0) {
		RELEASESHAREDB(res, array, f, cn);
		PyErr_SetString(PyExc_TypeError, "computeHessian: only for NGons.");
		return NULL;
	}

	posx = K_ARRAY::isCoordinateXPresent(varString);
	posy = K_ARRAY::isCoordinateYPresent(varString);
	posz = K_ARRAY::isCoordinateZPresent(varString);
	if (posx == -1 || posy == -1 || posz == -1) {
		PyErr_SetString(PyExc_TypeError, "computeHessian: coordinates not found in array.");
    	RELEASESHAREDB(res,array,f,cn);
		return NULL;
  	}
	
	posx++; posy++; posz++;

	// Check arrayc
  	char *varStringc; char *eltTypec;
  	FldArrayF *fc; FldArrayI *cnc;
  	E_Int nic, njc, nkc; // number of points of array
  	res = K_ARRAY::getFromArray(arrayc, varStringc, fc, nic, njc, nkc, cnc, eltTypec, true);

	E_Int npts = f->getSize();

	// Nombre de variables dont il faut calculer la hessienne
	E_Int nfld = fc->getNfld();
	vector<char *> vars;
	K_ARRAY::extractVars(varStringc, vars);

	vector<char *> varStrings;
	for (E_Int i = 0; i < nfld; i++) {
		char *local;
		computeHessianVarsString(vars[i], local);
		varStrings.push_back(local);
	}

	E_Int size = 0;
	for (E_Int i = 0; i < nfld; i++) {
		size += strlen(varStrings[i]) + 1;
	}

	char *varStringOut = new char [size];
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

	pt--; *pt = '\0';

	for (E_Int i = 0; i < vars.size(); i++) delete [] vars[i];

	// Calcul FE
	FldArrayI cFE;
	K_CONNECT::connectNG2FE(*cn, cFE);
	E_Int *cFE1 = cFE.begin(1);
	E_Int *cFE2 = cFE.begin(2);

	// Voisinage de chaque cellule
	FldArrayI cNFace;
	E_Int nelts;
	K_CONNECT::connectFE2NFace(cFE, cNFace, nelts);

	vector<vector<E_Int>> cell2cells(nelts);

	for (E_Int i = 0; i < nelts; i++)
		cell2cells[i].push_back(i);

	E_Int *cnp = cn->begin();
	E_Int nfaces = cnp[0];

	E_Int own, nei;
	for (E_Int i = 0; i < nfaces; i++) {
		own = cFE1[i]-1;
		nei = cFE2[i]-1;

		if (nei < 0) {
			continue;
		}
		
		cell2cells[own].push_back(nei);
		cell2cells[nei].push_back(own);
	}

	// Calcul des coordonnees des centres des faces := moyenne des coordonnes des noeuds
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

	// Calcul des coordonnees des centres des cellules := moyenne des coordonnes des centres des faces
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

	// Voisinage du voisinage des cellules
	E_Int nneis, NEI;
	cFE1 = cFE.begin(1);
	cFE2 = cFE.begin(2);
	std::unordered_set<E_Int> voisins;
	for (E_Int i = 0; i < nelts; i++) {
		nneis = cell2cells[i].size();
		voisins.clear();
		for (E_Int j = 1; j < nneis; j++) {
			cnp = cn->begin();
			nei = cell2cells[i][j];
			pos = posCells[nei];
			cnp += pos;
			stride = cnp[0];
			for (E_Int k = 1; k <= stride; k++) {
				face = cnp[k]-1;
				own = cFE1[face]-1;
				NEI = cFE2[face]-1;
				// bface: nth to do
				if (own < 0 || NEI < 0) continue;
				if (own == i || NEI == i) continue;
				if (own == nei)
					voisins.insert(NEI);
				else
					voisins.insert(own);
			}
		}
		for (auto it = voisins.begin(); it != voisins.end(); it++)
			cell2cells[i].push_back(*it);
	}

	// Construction de la matrice des moindres carres d'ordre 2
	E_Int nPoly = 6;
	E_Float xc, yc, zc, xci, yci, zci;
	E_Float dx, dy, dz;

	E_Float **beta = (E_Float **)malloc(nelts * sizeof(E_Float *));

	E_Float **A = (E_Float **)malloc(MAXNEIS * sizeof(E_Float *));
	for (E_Int i = 0; i < MAXNEIS; i++)
		A[i] = (E_Float *)malloc(nPoly * sizeof(E_Float));

	E_Float **tA = (E_Float **)malloc(nPoly * sizeof(E_Float *));
	for (E_Int i = 0; i < nPoly; i++)
		tA[i] = (E_Float *)malloc(MAXNEIS * sizeof(E_Float));

	E_Float *tAA = (E_Float *)malloc(nPoly * nPoly * sizeof(E_Float));

	E_Float *b = (E_Float *)malloc(MAXNEIS * sizeof(E_Float));

	E_Float *B = (E_Float *)malloc(nPoly * sizeof(E_Float));

	for (E_Int i = 0; i < nelts; i++) {	
		beta[i] = (E_Float *)malloc(nPoly * sizeof(E_Float));
		nneis = cell2cells[i].size();

		assert(nneis < MAXNEIS);
		
		xci = cellCenters(i, 1);
		yci = cellCenters(i, 2);
		zci = cellCenters(i, 3);

		// une equation par voisin
		for (E_Int j = 0; j < nneis; j++) {
			nei = cell2cells[i][j];
			xc = cellCenters(nei, 1);
			yc = cellCenters(nei, 2);
			zc = cellCenters(nei, 3);
			
			dx = xc-xci;
			dy = yc-yci;
			dz = zc-zci;

			A[j][0] = 1.;
			A[j][1] = dx;
			A[j][2] = dy;
			A[j][3] = dx*dx;
			A[j][4] = dx*dy;
			A[j][5] = dy*dy;
		}

		for (E_Int k = 0; k < nPoly; k++) {
			for (E_Int j = 0; j < nneis; j++) {
				tA[k][j] = A[j][k];
			}
		}	

		E_Int idx = 0;
		for (E_Int j = 0; j < nPoly; j++) {
			for (E_Int k = 0; k < nPoly; k++) {
				tAA[idx] = 0.;
				for (E_Int l = 0; l < nneis; l++) {
					tAA[idx] += tA[j][l] * A[l][k];
				}
				idx++;
			}
		}

		// b vector
		E_Float *fcp = fc->begin(1);
		for (E_Int j = 0; j < nneis; j++) {
			nei = cell2cells[i][j];
			b[j] = fcp[nei];
		}

		// B = tA*b
		for (E_Int j = 0; j < nPoly; j++) {
			B[j] = 0;
			for (E_Int k = 0; k < nneis; k++) {
				B[j] += tA[j][k] * b[k];
			}
		}

		// solve tAA.x = tA*b
		E_Int solved = K_LINEAR::solve(nPoly, 1, tAA, B, beta[i], "gauss"); 
	}

	// build
	PyObject *tpl = K_ARRAY::buildArray(9, varStringOut, npts,
		nelts, -1, eltType, true,
		cn->getSize()*cn->getNfld());

	E_Int *cnnp = K_ARRAY::getConnectPtr(tpl);

	K_KCORE::memcpy__(cnnp, cn->begin(), cn->getSize()*
		cn->getNfld());

	E_Float *Hp = K_ARRAY::getFieldPtr(tpl);

	FldArrayF H(nelts, 9, Hp, true);
	H.setAllValuesAtNull();

	E_Float *Hxx = H.begin(1);
	E_Float *Hxy = H.begin(2);
	E_Float *Hxz = H.begin(3);
	E_Float *Hyx = H.begin(4);
	E_Float *Hyy = H.begin(5);
	E_Float *Hyz = H.begin(6);
	E_Float *Hzx = H.begin(7);
	E_Float *Hzy = H.begin(8);
	E_Float *Hzz = H.begin(9);
	for (E_Int i = 0; i < nelts; i++) {
		Hxx[i] = 2.*beta[i][3];
		Hxy[i] = beta[i][4];
		Hxz[i] = 0.;
		Hyx[i] = beta[i][4];
		Hyy[i] = 2.*beta[i][5];
		Hyz[i] = 0.;
		Hzx[i] = 0.;
		Hzy[i] = 0.;
		Hzz[i] = 0.;
	}

	RELEASESHAREDB(res, array, f, cn);
	RELEASESHAREDB(res, arrayc, fc, cnc);

	delete [] varStringOut;
	for (E_Int i = 0; i < nelts; i++)
		free(beta[i]);
	free(beta);
	for (E_Int i = 0; i < MAXNEIS; i++)
		free(A[i]);
	free(A);
	for (E_Int i = 0; i < nPoly; i++)
		free(tA[i]);
	free(tA);
	free(tAA);
	free(b);
	free(B);

	return tpl;
}
