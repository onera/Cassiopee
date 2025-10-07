/*
    Copyright 2013-2025 Onera.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/

// application of a Gaussian Actuator Line -
// from Ronan Boisard Python source - RotorROM
# include "initiator.h"

using namespace K_FLD;
using namespace std;

#define RELEASEARRAYS                                       \
for (E_Int no = 0; no < nzones; no++)                       \
    RELEASESHAREDA(resl[no],objs[no],fields[no],a2[no],a3[no],a4[no]);

#define RELEASEROTMAT                                       \
    for (E_Int nob = 0; nob < NbBlades; nob++)              \
        for (E_Int nosec = 0; nosec < NbPoints; nosec++)    \
        {                                                   \
            E_Int ind = nosec + nob*NbPoints;               \
            PyObject* tpl = PyList_GetItem(listOfRotMat2SourceFrame, ind);\
            RELEASESHAREDN(tpl, vectOfRotMat2LocalFrame[ind]);\
        }

#define RELEASEBLOCK0                                       \
    for (E_Int j = 0; j < NbBlades; j++)                    \
    {                                                       \
        PyObject* tpl = PyList_GetItem(listOfAllLoads, j);  \
        RELEASESHAREDN(tpl, vectOfLoads[j]);                \
    }

#define RELEASEBLOCKP0                                      \
    for (E_Int j = 0; j < i; j++)                           \
    {                                                       \
        PyObject* tpl = PyList_GetItem(listOfAllLoads, j);  \
        RELEASESHAREDN(tpl, vectOfLoads[j]);                \
    }

#define RELEASEBLOCK1                                       \
    for (E_Int j = 0; j < NbBlades; j++)                    \
    {                                                       \
        PyObject* tpl = PyList_GetItem(listOfALPositions, j);         \
        RELEASESHAREDN(tpl, vectOfALPositions[j]);          \
    }

#define RELEASEBLOCKP1                                          \
    for (E_Int j = 0; j < i; j++)                               \
    {                                                           \
        PyObject* tpl = PyList_GetItem(listOfALPositions, j);   \
        RELEASESHAREDN(tpl, vectOfALPositions[j]);              \
    }

#define RELEASEBLOCK2                               \
    RELEASESHAREDN(pyLocalEpsX, localEps_X);  \
    RELEASESHAREDN(pyLocalEpsY, localEps_Y);  \
    RELEASESHAREDN(pyLocalEpsZ, localEps_Z);

// ============================================================================
/* AL Gaussian - version non normalisee (GaussInteg=1)*/
// ============================================================================
PyObject* K_INITIATOR::applyGaussianAL(PyObject* self, PyObject* args)
{
    char* TRUNCVARLOADS;
    char* TRUNCVARVELOS;
    E_Int NbBlades, NbPoints;
    PyObject *listOfAllLoads;// liste des numpy (3, nbsec) des composantes des forces pour chaque pale
    PyObject *pyLocalEpsX, *pyLocalEpsY, *pyLocalEpsZ;//local epsilon dans chaque direction pour les pts de la section
    PyObject *listOfALPositions;//coordonnees de tous les pts des actuator lines de toutes les pales
    PyObject* arrays; // champs en centres + coordonnees en centres
    PyObject* listOfRotMat2SourceFrame;// matrice de rotation pour chaque pale & pour chaque section de pale
    // la liste est de longueur nbblades*nbPts : toutes les sections de la 1ere pale, puis toutes les sections de la 2e etc
    if (!PYPARSETUPLE_(args, OOOO_ OOO_ II_ SS_,
                        &arrays, &listOfAllLoads, &listOfALPositions,
                        &listOfRotMat2SourceFrame,
                        &pyLocalEpsX, &pyLocalEpsY, &pyLocalEpsZ,
                        &NbBlades, &NbPoints, &TRUNCVARLOADS, &TRUNCVARVELOS))
        return NULL;

    vector<E_Int> resl;  vector<char*> varString;
    vector<FldArrayF*> fields;
    vector<void*> a2; //ni,nj,nk ou cnt en NS
    vector<void*> a3; //eltType en NS
    vector<void*> a4;
    vector<PyObject*> objs;
    E_Bool skipNoCoord = false; E_Bool skipStructured = false;
    E_Bool skipUnstructured = false; E_Bool skipDiffVars = false;
    E_Int ok = K_ARRAY::getFromArrays(arrays, resl, varString, fields, a2, a3, a4, objs,
                                      skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
    E_Int nzones = objs.size();
    if (ok == -1 || nzones == 0)
    {
        RELEASEARRAYS;
        PyErr_SetString(PyExc_TypeError,
                        "applyGaussianAL: 1st argument is not valid.");
        return NULL;
    }

    //recuperation des composantes des efforts sur les pales
    vector<FldArrayF*> vectOfLoads(NbBlades);
    for (E_Int i = 0; i < NbBlades; i++)
    {
        PyObject* pyLoadsForBlade = PyList_GetItem(listOfAllLoads, i);
        FldArrayF* loadsF;
        E_Int resn = K_NUMPY::getFromNumpyArray(pyLoadsForBlade, loadsF);
        if (resn == 0)
        {
            RELEASEARRAYS;
            RELEASEBLOCKP0;
            PyErr_SetString(PyExc_TypeError,
                    "applyAL: 2nd arg (loads per blade) must be a numpy of floats.");
            return NULL;
        }
        vectOfLoads[i] = loadsF;
    }
    //recuperation des positions des AL de toutes les pales
    vector<FldArrayF*> vectOfALPositions(NbBlades);
    for (E_Int i = 0; i < NbBlades; i++)
    {
        PyObject* pyALPositionsForBlade = PyList_GetItem(listOfALPositions, i);
        FldArrayF* ALPositionsF;
        E_Int resn = K_NUMPY::getFromNumpyArray(pyALPositionsForBlade, ALPositionsF);
        if (resn == 0)
        {
            RELEASEARRAYS;
            RELEASEBLOCK0;
            RELEASEBLOCKP1;
            PyErr_SetString(PyExc_TypeError,
                            "applyAL: 3rd arg (AL positions per blade) must be a numpy of floats.");
            return NULL;
        }
        vectOfALPositions[i] = ALPositionsF;
    }
    vector< FldArrayF*> vectOfRotMat2LocalFrame(NbBlades*NbPoints);
    for (E_Int nob = 0; nob < NbBlades; nob++)
        for (E_Int nosec = 0; nosec < NbPoints; nosec++)
        {
            E_Int ind = nosec + nob*NbPoints;
            PyObject* tpl = PyList_GetItem(listOfRotMat2SourceFrame, ind);
            FldArrayF* RotMatLocal;
            E_Int resn = K_NUMPY::getFromNumpyArray(tpl, RotMatLocal);
            if (resn==0)
            {
                RELEASEARRAYS;
                RELEASEBLOCK0;
                RELEASEBLOCK1;
                for (E_Int ind2 = 0; ind2 < ind; ind2++)
                {
                    PyObject* tpl = PyList_GetItem(listOfRotMat2SourceFrame, ind2);
                    RELEASESHAREDN(tpl, vectOfRotMat2LocalFrame[ind]);
                }
            }
            vectOfRotMat2LocalFrame[ind] = RotMatLocal;
        }

    // recuperation des epsilons locaux a chq section de pale
    // chaque pale a le meme epsx,epsy,epsz pour une section donnee
    FldArrayF* localEps_X;
    E_Int rese = K_NUMPY::getFromNumpyArray(pyLocalEpsX, localEps_X);
    if (rese == 0)
    {
        RELEASEARRAYS;
        RELEASEBLOCK0;
        RELEASEBLOCK1; RELEASEROTMAT;
        PyErr_SetString(PyExc_TypeError,
                        "applyAL: 3rd arg (eps local wrt x) must be a numpy of floats .");
        return NULL;
    }

    FldArrayF* localEps_Y;
    rese = K_NUMPY::getFromNumpyArray(pyLocalEpsY, localEps_Y);
    if (rese == 0)
    {
        RELEASEARRAYS;
        RELEASEBLOCK0;
        RELEASEBLOCK1; RELEASEROTMAT;
        RELEASESHAREDN(pyLocalEpsX, localEps_X);
        PyErr_SetString(PyExc_TypeError,
                        "applyAL: 4th arg (eps local wrt y) must be a numpy of floats .");
        return NULL;
    }

    FldArrayF* localEps_Z;
    rese = K_NUMPY::getFromNumpyArray(pyLocalEpsZ, localEps_Z);
    if (rese == 0)
    {
        RELEASEARRAYS;
        RELEASEBLOCK0;
        RELEASEBLOCK1; RELEASEROTMAT;
        RELEASESHAREDN(pyLocalEpsX, localEps_X);
        RELEASESHAREDN(pyLocalEpsY, localEps_Y);
        PyErr_SetString(PyExc_TypeError,
                        "applyAL: 5th arg (eps local wrt z) must be a numpy of floats .");
        return NULL;
    }

    E_Float* epsX = localEps_X->begin();
    E_Float* epsY = localEps_Y->begin();
    E_Float* epsZ = localEps_Z->begin();
    FldArrayF epsXInvF(NbPoints); E_Float* epsXInv = epsXInvF.begin();
    FldArrayF epsYInvF(NbPoints); E_Float* epsYInv = epsYInvF.begin();
    FldArrayF epsZInvF(NbPoints); E_Float* epsZInv = epsZInvF.begin();
    FldArrayF eps3InvF(NbPoints);  E_Float* eps3Inv = eps3InvF.begin();
 #pragma omp parallel default(shared)
  {
#pragma omp for
    for (E_Int nosec = 0; nosec < NbPoints; nosec++)
    {
        epsXInv[nosec]= 1./epsX[nosec];
        epsYInv[nosec]= 1./epsY[nosec];
        epsZInv[nosec]= 1./epsZ[nosec];
        eps3Inv[nosec]= epsXInv[nosec]*epsYInv[nosec]*epsZInv[nosec];
    }
  }
    E_Float GaussianNormInv = 1.;//on ne normalise pas la gaussienne - cas le plus simple
    E_Float alphaPi = 1./sqrt(K_CONST::E_PI*K_CONST::E_PI*K_CONST::E_PI);

    for (E_Int nob = 0; nob < NbBlades; nob++)
    {
        E_Float* x_AL = vectOfALPositions[nob]->begin(1);
        E_Float* y_AL = vectOfALPositions[nob]->begin(2);
        E_Float* z_AL = vectOfALPositions[nob]->begin(3);
        E_Float* Fx = vectOfLoads[nob]->begin();
        E_Float* Fy = vectOfLoads[nob]->begin()+NbPoints;
        E_Float* Fz = vectOfLoads[nob]->begin()+2*NbPoints;

        for(E_Int nosec = 0; nosec < NbPoints; nosec++)
        {
            E_Float Fx0 = Fx[nosec];
            E_Float Fy0 = Fy[nosec];
            E_Float Fz0 = Fz[nosec];
            E_Float epsxi = epsXInv[nosec];
            E_Float epsyi = epsYInv[nosec];
            E_Float epszi = epsZInv[nosec];

            E_Int nomat = nosec + nob * NbPoints;
            FldArrayF& RotMatLocal = *(vectOfRotMat2LocalFrame[nomat]);
            E_Float RotMat11 = epsxi*RotMatLocal(0,1);
            E_Float RotMat12 = epsyi*RotMatLocal(1,1);
            E_Float RotMat13 = epszi*RotMatLocal(2,1);
            E_Float RotMat21 = epsxi*RotMatLocal(0,2);
            E_Float RotMat22 = epsyi*RotMatLocal(1,2);
            E_Float RotMat23 = epszi*RotMatLocal(2,2);
            E_Float RotMat31 = epsxi*RotMatLocal(0,3);
            E_Float RotMat32 = epsyi*RotMatLocal(1,3);
            E_Float RotMat33 = epszi*RotMatLocal(2,3);

            E_Float coef = alphaPi * eps3Inv[nosec];

            for (E_Int noz = 0; noz < nzones; noz++)
            {
                E_Int posx = K_ARRAY::isNamePresent("XCell", varString[noz]); //K_ARRAY::isCoordinateXPresent(varString);
                E_Int posy = K_ARRAY::isNamePresent("YCell", varString[noz]);//K_ARRAY::isCoordinateYPresent(varString);
                E_Int posz = K_ARRAY::isNamePresent("ZCell", varString[noz]);//K_ARRAY::isCoordinateZPresent(varString);
                E_Int posvol = K_ARRAY::isNamePresent("vol", varString[noz]);
                E_Int postruncl = K_ARRAY::isNamePresent(TRUNCVARLOADS, varString[noz]);
                E_Int postruncv = K_ARRAY::isNamePresent(TRUNCVARVELOS, varString[noz]);
                E_Int posgs = K_ARRAY::isNamePresent("GaussSum", varString[noz]);
                E_Int posrou_src = K_ARRAY::isNamePresent("MomentumX_src", varString[noz]);
                E_Int posrov_src = K_ARRAY::isNamePresent("MomentumY_src", varString[noz]);
                E_Int posrow_src = K_ARRAY::isNamePresent("MomentumZ_src", varString[noz]);
                E_Int posroe_src = K_ARRAY::isNamePresent("EnergyStagnationDensity_src", varString[noz]);
                E_Int posu = K_ARRAY::isVelocityXPresent(varString[noz]);
                E_Int posv = K_ARRAY::isVelocityYPresent(varString[noz]);
                E_Int posw = K_ARRAY::isVelocityZPresent(varString[noz]);
                E_Int posuo = K_ARRAY::isNamePresent("VelXOut", varString[noz]);
                E_Int posvo = K_ARRAY::isNamePresent("VelYOut", varString[noz]);
                E_Int poswo = K_ARRAY::isNamePresent("VelZOut", varString[noz]);
                if (posx==-1 || posy==-1 || posz==-1 || posvol==-1 ||
                    postruncl==-1 || postruncv==-1 || posgs==-1 ||
                    posrou_src==-1 || posrov_src==-1 || posrow_src==-1 || posroe_src == -1 ||
                    posu==-1 || posv==-1 || posw==-1 ||
                    posuo==-1 || posvo==-1 || poswo==-1)
                {
                    RELEASEARRAYS;
                    RELEASEBLOCK0;
                    RELEASEBLOCK1; RELEASEROTMAT;
                    RELEASESHAREDN(pyLocalEpsX, localEps_X);
                    RELEASESHAREDN(pyLocalEpsY, localEps_Y);
                    RELEASESHAREDN(pyLocalEpsZ, localEps_Z);

                    printf(" posx = " SF_D3_ " | posvol " SF_D_ " postruncl " SF_D_ " postruncv " SF_D_ " | posgs = " SF_D_ " posrou_src  = (" SF_D4_ ") | posu = " SF_D3_ " | posuo = " SF_D3_ " \n",
                            posx, posy, posz, posvol, postruncl, postruncv, posgs,
                            posrou_src, posrov_src, posrow_src, posroe_src, posu, posv, posw, posuo, posvo, poswo );
                    PyErr_SetString(PyExc_TypeError,
                    "applyGaussianAL: cannot find required fields in array.");
                    return NULL;
                }
                posx++; posy++; posz++; posvol++;
                postruncl++; postruncv++; posgs++;
                posrou_src++; posrov_src++; posrow_src++;
                posuo++; posvo++; poswo++;
                FldArrayF* f = fields[noz];
                E_Int npts = f->getSize();
                E_Float* xt = f->begin(posx);
                E_Float* yt = f->begin(posy);
                E_Float* zt = f->begin(posz);
                E_Float* volt = f->begin(posvol);
                E_Float* truncloadst = f->begin(postruncl);
                E_Float* truncvelost = f->begin(postruncv);
                E_Float* gausssumt = f->begin(posgs);
                E_Float* MomentumX_src = f->begin(posrou_src);
                E_Float* MomentumY_src = f->begin(posrov_src);
                E_Float* MomentumZ_src = f->begin(posrow_src);
                E_Float* EnergyStagnationDensity_src = f->begin(posroe_src);
                E_Float* vxt = f->begin(posu);
                E_Float* vyt = f->begin(posv);
                E_Float* vzt = f->begin(posw);
                E_Float* vxo = f->begin(posuo);
                E_Float* vyo = f->begin(posvo);
                E_Float* vzo = f->begin(poswo);

            #pragma omp parallel default(shared)
                {
                #pragma omp for
                for (E_Int ind = 0; ind < npts; ind++)
                {
                    E_Float volCell = volt[ind];
                    E_Float truncLoadsCell = truncloadst[ind];
                    E_Float truncVelosCell = truncvelost[ind];
                    E_Float XCell = xt[ind];
                    E_Float YCell = yt[ind];
                    E_Float ZCell = zt[ind];

                    //(Xcell-X_AL)/epsX
                    E_Float X_ind =
                        RotMat11*(XCell-x_AL[nosec])+
                        RotMat12*(YCell-y_AL[nosec])+
                        RotMat13*(ZCell-z_AL[nosec]);

                    E_Float Y_ind =
                        RotMat21*(XCell-x_AL[nosec])+
                        RotMat22*(YCell-y_AL[nosec])+
                        RotMat23*(ZCell-z_AL[nosec]);

                    E_Float Z_ind =
                        RotMat31*(XCell-x_AL[nosec])+
                        RotMat32*(YCell-y_AL[nosec])+
                        RotMat33*(ZCell-z_AL[nosec]);

                    E_Float gauss = coef * volCell * exp(-X_ind*X_ind-Y_ind*Y_ind-Z_ind*Z_ind);
                    E_Float gaussl = gauss * truncLoadsCell * GaussianNormInv;
                    E_Float gaussFx = -gaussl*Fx0;
                    E_Float gaussFy = -gaussl*Fy0;
                    E_Float gaussFz = -gaussl*Fz0;
                    MomentumX_src[ind] += gaussFx;
                    MomentumY_src[ind] += gaussFy;
                    MomentumZ_src[ind] += gaussFz;
                    EnergyStagnationDensity_src[ind]+=gaussFx*vxt[ind]+gaussFy*vyt[ind]+gaussFz*vzt[ind];
                    gausssumt[ind]+=gaussl;
                    E_Float gaussvel = gauss * truncVelosCell * GaussianNormInv;
                    vxo[ind] = vxt[ind]*gaussvel;
                    vyo[ind] = vyt[ind]*gaussvel;
                    vzo[ind] = vzt[ind]*gaussvel;
                }//npts
                }//omp
            }//nzones
        }//nsections
    }//nblades
// for (E_Int noz = 0; noz < nzones; noz++)
// {
//     FldArrayF* f = fields[noz];
//     E_Int npts = f->getSize();
//     E_Int posrou_src = K_ARRAY::isNamePresent("MomentumX_src", varString[noz]);
//     E_Int posrov_src = K_ARRAY::isNamePresent("MomentumY_src", varString[noz]);
//     E_Int posrow_src = K_ARRAY::isNamePresent("MomentumZ_src", varString[noz]);
//     E_Int posu = K_ARRAY::isVelocityXPresent(varString[noz]);
//     E_Int posv = K_ARRAY::isVelocityYPresent(varString[noz]);
//     E_Int posw = K_ARRAY::isVelocityZPresent(varString[noz]);
//     E_Int posroe_src = K_ARRAY::isNamePresent("EnergyStagnationDensity_src", varString[noz]);
//     posu++; posv++; posw++;
//     posrou_src++; posrov_src++; posrow_src++; posroe_src++;
//     E_Float* vxt = f->begin(posu);
//     E_Float* vyt = f->begin(posv);
//     E_Float* vzt = f->begin(posw);
//     E_Float* MomentumX_src = f->begin(posrou_src);
//     E_Float* MomentumY_src = f->begin(posrov_src);
//     E_Float* MomentumZ_src = f->begin(posrow_src);
//     E_Float* EnergyStagnationDensity_src = f->begin(posroe_src);

// #pragma omp parallel default(shared)
// {
// #pragma omp for
//     for (E_Int ind = 0; ind < npts; ind++)
//     {
//         EnergyStagnationDensity_src[ind]=MomentumX_src[ind]*vxt[ind]+MomentumY_src[ind]*vyt[ind]+MomentumZ_src[ind]*vzt[ind];
//     }
// }
// }
    RELEASEARRAYS;
    RELEASEBLOCK0;
    RELEASEBLOCK1;
    RELEASEBLOCK2;
    RELEASEROTMAT;
    Py_INCREF(Py_None);
    return Py_None;
}

