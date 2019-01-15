/*    
    Copyright 2013-2019 Onera.

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

# include "connector.h"
# include "Search/BbTree.h"
# include "Search/KdTree.h"
# include "Fld/ArrayAccessor.h"
using namespace K_FLD;
using namespace std;
using namespace K_SEARCH;

#define RELEASEZONES \
    for (E_Int nos = 0; nos < nzonesS; nos++)\
        RELEASESHAREDS(objst[nos], structF[nos]);\
    for (E_Int nos = 0; nos < nzones; nos++)\
        RELEASESHAREDU(objut[nos], unstrF[nos], cnt[nos]);

#define RELEASEBODIES \
    for (E_Int nos = 0; nos < nbodiesS; nos++)\
        RELEASESHAREDS(objsb[nos], structbF[nos]);\
    for (E_Int nos = 0; nos < nbodies; nos++)\
        RELEASESHAREDU(objub[nos], unstrbF[nos], cnb[nos]);

#define RELEASEFRONT\
    for (E_Int nos = 0; nos < nfrontS; nos++)\
        RELEASESHAREDS(objsf[nos], structfF[nos]);\
    for (E_Int nos = 0; nos < nfronts; nos++)\
        RELEASESHAREDU(objuf[nos], unstrfF[nos], cnf[nos]);


// ============================================================================
    /* Get the wall and image points if front is provided 
    wall pts are obtained by projDir (projOrtho, closest pt if impossible)
    image pts are also obtained by projections*/
// ============================================================================
PyObject* K_CONNECTOR::getIBMPtsWithFront(PyObject* self, PyObject* args)
{
    PyObject *allCorrectedPts, *bodySurfaces, *frontSurfaces, *normalNames;
    PyObject *ListOfSnearsLoc;
    E_Int signOfDist; //if correctedPts are inside bodies: sign = -1, else sign=1
    E_Int depth;//nb of layers of ghost cells
    if (!PYPARSETUPLEI(args,"OOOOOll","OOOOOii", 
                       &allCorrectedPts, &ListOfSnearsLoc, &bodySurfaces, &frontSurfaces, 
                       &normalNames, &signOfDist, &depth)) return NULL;

    // extract list of snearsloc
    if (PyList_Size(ListOfSnearsLoc) == 0)
    {
        PyErr_SetString(PyExc_TypeError, 
                        "getIBMPtsWithFront: 2nd argument is an empty list.");
        return NULL;
    }


    E_Int nsnears = PyList_Size(ListOfSnearsLoc);
    vector<E_Float> vectOfSnearsLoc;
    PyObject* tpl0 = NULL;

    for (int i = 0; i < nsnears; i++)
    {
        tpl0 = PyList_GetItem(ListOfSnearsLoc, i);
        if (PyFloat_Check(tpl0) == 0)
        {
            PyErr_SetString(PyExc_TypeError, 
                            "getIBMPtsWithFront: not a valid value for snear.");
            return NULL;
        } 
        else vectOfSnearsLoc.push_back(PyFloat_AsDouble(tpl0));
    }
  
    E_Int sign = -signOfDist; // sens de projection sur la paroi

    //tolerance of distance between corrected pts and wall/image to control projection
    //distance max for image pt to its corrected pt : (depth+1)*sqrt(2)*snearloc 
    E_Float toldistFactorImage;
    if (signOfDist==-1)// Euler : we move away the front from one additional layer
        toldistFactorImage = (depth+2)*(depth+2)*2;
    else 
        toldistFactorImage = (depth+1)*(depth+1)*2;
    //distance max for wall pt to its corrected pt ( depth)*sqrt(2)*snearloc
    E_Float toldistFactorWall = (depth+1)*(depth+1)*2;

    // Check normal components
    if (PyList_Check(normalNames) == 0)
    {
        PyErr_SetString(PyExc_TypeError, 
                        "getIBMPts: normal vars must be a list of strings.");
        return NULL;
    }
    E_Int nvars = PyList_Size(normalNames);
    if (nvars != 3)
    {
        PyErr_SetString(PyExc_TypeError, 
                        "getIBMPts: normal vars must be a 3-component vector.");
        return NULL;
    }
    char* var;
    vector<char*> varsn;// normal components names
    for (E_Int v = 0; v < 3; v++)
    {
        if (PyString_Check(PyList_GetItem(normalNames, v)) == 0)
        {
            PyErr_SetString(PyExc_TypeError,
                            "getIBMPts: invalid string for normal component.");
            return NULL;
        }
        else 
        {
            var = PyString_AsString(PyList_GetItem(normalNames, v));
            varsn.push_back(var);
        }
    }
    // Extract correctedPts
    vector<E_Int> resl;
    vector<char*> structVarString; vector<char*> unstrVarString;
    vector<FldArrayF*> structF; vector<FldArrayF*> unstrF;
    vector<E_Int> nit; vector<E_Int> njt; vector<E_Int> nkt;
    vector<FldArrayI*> cnt; vector<char*> eltType;
    vector<PyObject*> objst, objut;
    E_Boolean skipNoCoord = true;
    E_Boolean skipStructured = true;
    E_Boolean skipUnstructured = false;
    E_Boolean skipDiffVars = false;
    E_Int isOk = K_ARRAY::getFromArrays(allCorrectedPts, resl, structVarString, unstrVarString,
                                        structF, unstrF, nit, njt, nkt, cnt, eltType, objst, objut, 
                                        skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
    E_Int nzones = objut.size(); E_Int nzonesS = objst.size();
    if (isOk == -1)   
    {
        PyErr_SetString(PyExc_TypeError,"getIBMPts: 1st arg is not valid.");
        RELEASEZONES;
        return NULL;
    }
    if ( nzones != nsnears)
    {
        PyErr_SetString(PyExc_TypeError,"getIBMPts: 1st and 2nd arg must be lists of same length.");
        RELEASEZONES; return NULL;      
    }
    for (E_Int no = 0; no < nzones; no++)
    {
        if (K_STRING::cmp(eltType[no], "NODE") != 0) 
        {
            PyErr_SetString(PyExc_TypeError,"getIBMPts: 1st arg must be a NODE type zone.");
            RELEASEZONES;
            return NULL;
        }
    }
    // corrected points coordinates and normal components: position 
    E_Int posx1, posy1, posz1;
    vector<E_Int> posxt; vector<E_Int> posyt; vector<E_Int> poszt;
    vector<E_Int> posnxt; vector<E_Int> posnyt; vector<E_Int> posnzt;
    for (E_Int no = 0; no < nzones; no++)        
    {
        posx1 = K_ARRAY::isCoordinateXPresent(unstrVarString[no]); posx1++;
        posy1 = K_ARRAY::isCoordinateYPresent(unstrVarString[no]); posy1++;
        posz1 = K_ARRAY::isCoordinateZPresent(unstrVarString[no]); posz1++;
        posxt.push_back(posx1); posyt.push_back(posy1); poszt.push_back(posz1);

        posx1 = K_ARRAY::isNamePresent(varsn[0], unstrVarString[no]);      
        posy1 = K_ARRAY::isNamePresent(varsn[1], unstrVarString[no]);      
        posz1 = K_ARRAY::isNamePresent(varsn[2], unstrVarString[no]);           
        if (posx1==-1 || posy1==-1 || posz1==-1) 
        {
            PyErr_SetString(PyExc_TypeError,"getIBMPts: one of normal components not found in 1st arg.");
            RELEASEZONES; 
            return NULL;
        }
        else
        {
            posx1++; posy1++; posz1++;
            posnxt.push_back(posx1); posnyt.push_back(posy1); posnzt.push_back(posz1);
        }
    }

    // extract body surfaces
    vector<E_Int> resb;
    vector<char*> structVarStringb; vector<char*> unstrVarStringb;
    vector<FldArrayF*> structbF; vector<FldArrayF*> unstrbF;
    vector<E_Int> nib; vector<E_Int> njb; vector<E_Int> nkb;
    vector<FldArrayI*> cnb; vector<char*> eltTypeb;
    vector<PyObject*> objsb, objub;
    skipNoCoord = true;
    skipStructured = true;
    skipUnstructured = false;
    skipDiffVars = false;
    isOk = K_ARRAY::getFromArrays(bodySurfaces, resb, structVarStringb, unstrVarStringb,
                                  structbF, unstrbF, nib, njb, nkb, cnb, eltTypeb, objsb, objub, 
                                  skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
    E_Int nbodies = objub.size(); E_Int nbodiesS = objsb.size();
    if (isOk == -1)   
    {
        PyErr_SetString(PyExc_TypeError,"getIBMPts: 3rd arg is not valid.");
        RELEASEZONES; RELEASEBODIES;
        return NULL;
    }
    for (E_Int no = 0; no < nbodies; no++)
    {
        if (K_STRING::cmp(eltTypeb[no], "TRI") != 0) 
        {
            PyErr_SetString(PyExc_TypeError,"getIBMPts: 3rd arg must be a TRI type zone.");
            RELEASEZONES; RELEASEBODIES;
            return NULL;
        }
    }

    vector<E_Int> posxb; vector<E_Int> posyb; vector<E_Int> poszb;
    for (E_Int no = 0; no < nbodies; no++)
    {
        posx1 = K_ARRAY::isCoordinateXPresent(unstrVarStringb[no]); posx1++;
        posy1 = K_ARRAY::isCoordinateYPresent(unstrVarStringb[no]); posy1++;
        posz1 = K_ARRAY::isCoordinateZPresent(unstrVarStringb[no]); posz1++;
        posxb.push_back(posx1); posyb.push_back(posy1); poszb.push_back(posz1); 
    }

    // extract front surfaces
    vector<E_Int> resf;
    vector<char*> structVarStringf; vector<char*> unstrVarStringf;
    vector<FldArrayF*> structfF; vector<FldArrayF*> unstrfF;
    vector<E_Int> nif; vector<E_Int> njf; vector<E_Int> nkf;
    vector<FldArrayI*> cnf; vector<char*> eltTypef;
    vector<PyObject*> objsf, objuf;
    skipNoCoord = true;
    skipStructured = true;
    skipUnstructured = false;
    skipDiffVars = false;
    isOk = K_ARRAY::getFromArrays(frontSurfaces, resf, structVarStringf, unstrVarStringf,
                                  structfF, unstrfF, nif, njf, nkf, cnf, eltTypef, objsf, objuf, 
                                  skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
    E_Int nfronts = objuf.size(); E_Int nfrontS = objsf.size();
    if (isOk == -1)   
    {
        PyErr_SetString(PyExc_TypeError,"getIBMPts: 3rd arg is not valid.");
        RELEASEZONES; RELEASEBODIES; RELEASEFRONT;
        return NULL;
    }    
    for (E_Int no = 0; no < nfronts; no++)
    {
        if (K_STRING::cmp(eltTypef[no], "TRI") != 0) 
        {
            PyErr_SetString(PyExc_TypeError,"getIBMPts: 3rd arg must be a TRI type zone.");
            RELEASEZONES; RELEASEBODIES; RELEASEFRONT;
            return NULL;
        }
    }

    vector<E_Int> posxf; vector<E_Int> posyf; vector<E_Int> poszf;
    for (E_Int no = 0; no < nfronts; no++)
    {
        posx1 = K_ARRAY::isCoordinateXPresent(unstrVarStringf[no]); posx1++;
        posy1 = K_ARRAY::isCoordinateYPresent(unstrVarStringf[no]); posy1++;
        posz1 = K_ARRAY::isCoordinateZPresent(unstrVarStringf[no]); posz1++;
        posxf.push_back(posx1); posyf.push_back(posy1); poszf.push_back(posz1); 
    }

    // Checks completed...build arrays
    PyObject* PyListOfImagePts = PyList_New(0);  
    PyObject* PyListOfWallPts  = PyList_New(0);
    
    E_Int nvarOut = 3;
    char varStringOut[K_ARRAY::VARSTRINGLENGTH]; 
    strcpy(varStringOut,"CoordinateX,CoordinateY,CoordinateZ");
    char eltTypeOut[8]; strcpy(eltTypeOut,"NODE");
    FldArrayI cnl(0);
    E_Int nelts = 0;
    vector<E_Float*> xit; vector<E_Float*> yit; vector<E_Float*> zit;
    vector<E_Float*> xwt; vector<E_Float*> ywt; vector<E_Float*> zwt;

    for (E_Int no = 0; no < nzones; no++)
    {
        E_Int npts = unstrF[no]->getSize();
        PyObject* tpli = K_ARRAY::buildArray(nvarOut, varStringOut, npts, nelts, -1, eltTypeOut, false);
        E_Float* coordip = K_ARRAY::getFieldPtr(tpli);
        FldArrayF coordi(npts, 3, coordip, true); // non initialized yet
        coordi.setOneField(*unstrF[no],posxt[no],1);
        coordi.setOneField(*unstrF[no],posyt[no],2);
        coordi.setOneField(*unstrF[no],poszt[no],3);
        xit.push_back(coordi.begin(1));
        yit.push_back(coordi.begin(2));
        zit.push_back(coordi.begin(3));
        PyList_Append(PyListOfImagePts, tpli); Py_DECREF(tpli);

        PyObject* tplw = K_ARRAY::buildArray(nvarOut, varStringOut, npts, nelts, -1, eltTypeOut, false);
        E_Float* coordwp = K_ARRAY::getFieldPtr(tplw);
        FldArrayF coordw(npts, 3, coordwp, true); // non initialized yet
        coordw.setOneField(*unstrF[no],posxt[no],1);
        coordw.setOneField(*unstrF[no],posyt[no],2);
        coordw.setOneField(*unstrF[no],poszt[no],3);
        xwt.push_back(coordw.begin(1));
        ywt.push_back(coordw.begin(2));
        zwt.push_back(coordw.begin(3));
        PyList_Append(PyListOfWallPts, tplw); Py_DECREF(tplw);        
    }

    //Creation of BBTrees and kdtree for body surfaces
    typedef K_SEARCH::BoundingBox<3>  BBox3DType; 
    vector<K_SEARCH::BbTree3D*> vectOfBodyBBTrees;
    E_Float minB[3];  E_Float maxB[3];
    vector< vector<BBox3DType*> > vectOfBodyBoxes;// to be deleted at the end
    E_Int nptsBodiesMax=0;
    for (E_Int v = 0; v < nbodies; v++) nptsBodiesMax += unstrbF[v]->getSize();

    FldArrayF* bodySurfacesPts = new FldArrayF(nptsBodiesMax,3);
    E_Float* xb2 = bodySurfacesPts->begin(1);
    E_Float* yb2 = bodySurfacesPts->begin(2);
    E_Float* zb2 = bodySurfacesPts->begin(3);
    E_Int nopt=0;
    for (E_Int v = 0; v < nbodies; v++)
    {
        E_Float* xb = unstrbF[v]->begin(posxb[v]);
        E_Float* yb = unstrbF[v]->begin(posyb[v]);
        E_Float* zb = unstrbF[v]->begin(poszb[v]);

        for (E_Int ind = 0; ind < unstrbF[v]->getSize(); ind++)
        {
          xb2[nopt] = xb[ind]; yb2[nopt] = yb[ind]; zb2[nopt] = zb[ind];
          nopt++;
        }
        E_Int neltsb = cnb[v]->getSize();
        vector<BBox3DType*> boxes(neltsb);// list of bboxes of all triangles
        FldArrayF bbox(neltsb,6);// xmin, ymin, zmin, xmax, ymax, zmax 
        K_COMPGEOM::boundingBoxOfUnstrCells(*cnb[v], xb, yb, zb, bbox);
        E_Float* xminp = bbox.begin(1); E_Float* xmaxp = bbox.begin(4);
        E_Float* yminp = bbox.begin(2); E_Float* ymaxp = bbox.begin(5);
        E_Float* zminp = bbox.begin(3); E_Float* zmaxp = bbox.begin(6);
        for (E_Int et = 0; et < neltsb; et++)
        {
          minB[0] = xminp[et]; minB[1] = yminp[et]; minB[2] = zminp[et];
          maxB[0] = xmaxp[et]; maxB[1] = ymaxp[et]; maxB[2] = zmaxp[et]; 
          boxes[et] = new BBox3DType(minB, maxB);
        }
        vectOfBodyBoxes.push_back(boxes);

        // Build the box tree.
        K_SEARCH::BbTree3D* bbtree = new K_SEARCH::BbTree3D(boxes);
        vectOfBodyBBTrees.push_back(bbtree);
    }

    ArrayAccessor<FldArrayF> coordAcc(*bodySurfacesPts, 1,2,3);
    KdTree<FldArrayF> kdt(coordAcc, E_EPSILON);
    // Creation of BBTrees for front surfaces 
    vector<K_SEARCH::BbTree3D*> vectOfFrontBBTrees;
    vector< vector<BBox3DType*> > vectOfFrontBoxes;// to be deleted at the end
    E_Int nptsFrontMax=0;
    for (E_Int v = 0; v < nfronts; v++)
        nptsFrontMax += unstrfF[v]->getSize();
    FldArrayF* frontSurfacesPts = new FldArrayF(nptsFrontMax,3);
    E_Float* xf2 = frontSurfacesPts->begin(1);
    E_Float* yf2 = frontSurfacesPts->begin(2);
    E_Float* zf2 = frontSurfacesPts->begin(3);
    E_Int noptf=0;
    for (E_Int v = 0; v < nfronts; v++)
    {
        E_Float* xf = unstrfF[v]->begin(posxf[v]);
        E_Float* yf = unstrfF[v]->begin(posyf[v]);
        E_Float* zf = unstrfF[v]->begin(poszf[v]);
        for (E_Int ind = 0; ind < unstrfF[v]->getSize(); ind++)
        {
          xf2[noptf] = xf[ind]; yf2[noptf] = yf[ind]; zf2[noptf] = zf[ind];
          noptf++;
        }
        E_Int neltsf = cnf[v]->getSize();
        vector<BBox3DType*> boxes(neltsf);// list of bboxes of all triangles
        FldArrayF bbox(neltsf,6);// xmin, ymin, zmin, xmax, ymax, zmax 
        K_COMPGEOM::boundingBoxOfUnstrCells(*cnf[v], xf, yf, zf, bbox);

        E_Float* xminp = bbox.begin(1); E_Float* xmaxp = bbox.begin(4);
        E_Float* yminp = bbox.begin(2); E_Float* ymaxp = bbox.begin(5);
        E_Float* zminp = bbox.begin(3); E_Float* zmaxp = bbox.begin(6);
        for (E_Int et = 0; et < neltsf; et++)
        {
          minB[0] = xminp[et]; minB[1] = yminp[et]; minB[2] = zminp[et];
          maxB[0] = xmaxp[et]; maxB[1] = ymaxp[et]; maxB[2] = zmaxp[et]; 
          boxes[et] = new BBox3DType(minB, maxB);
        }
        vectOfFrontBoxes.push_back(boxes);

        // Build the box tree.
        K_SEARCH::BbTree3D* bbtree = new K_SEARCH::BbTree3D(boxes);
        vectOfFrontBBTrees.push_back(bbtree);
    }

    ArrayAccessor<FldArrayF> coordAccF(*frontSurfacesPts, 1,2,3);
    KdTree<FldArrayF> kdtf(coordAccF, E_EPSILON);
    /*--------------------------------------------------------*/
    // projectDir of all the points onto the bodies-> wall pts
    // projectDir of all the points onto the front -> image pts
    /*--------------------------------------------------------*/    
    E_Float tol = K_CONST::E_GEOM_CUTOFF;
    E_Float dirx0, diry0, dirz0, xsav, ysav, zsav;
    E_Float xc0, yc0, zc0, xw0, yw0, zw0, xi0, yi0, zi0;
    E_Float dist2, distl;
    vector<E_Int> indicesBB; 
    E_Float pr1[3]; E_Float pr2[3]; E_Float pt[3];
    E_Int oriented = 1;
    E_Int ok, notri, indvert1, indvert2, indvert3, indp, err;
    E_Float rx, ry, rz, rad;

    //distance of corrected pts to wall pts and image pts 
    E_Float delta1, delta2;
    //Corresponding directions
    //E_Float nx, ny, nz;
    E_Float nxp, nyp, nzp, nxs, nys, nzs;
    E_Int *cnVert1, *cnVert2, *cnVert3;
    E_Float edgeLen1, edgeLen2;
    E_Int foundw, foundi, valid;
    PyObject* PyListIndicesByIBCType = PyList_New(0);

    for (E_Int noz = 0; noz < nzones; noz++)
    {   
        FldArrayF* correctedPts = unstrF[noz];
        E_Int npts = correctedPts->getSize();
        E_Float* ptrNX = correctedPts->begin(posnxt[noz]);
        E_Float* ptrNY = correctedPts->begin(posnyt[noz]);
        E_Float* ptrNZ = correctedPts->begin(posnzt[noz]);
        E_Float* ptrXC = correctedPts->begin(posxt[noz]);
        E_Float* ptrYC = correctedPts->begin(posyt[noz]);
        E_Float* ptrZC = correctedPts->begin(poszt[noz]);
        E_Float* ptrXW = xwt[noz];
        E_Float* ptrYW = ywt[noz];
        E_Float* ptrZW = zwt[noz];
        E_Float* ptrXI = xit[noz];
        E_Float* ptrYI = yit[noz];
        E_Float* ptrZI = zit[noz];
        E_Float snearloc = vectOfSnearsLoc[noz];        

        vector<FldArrayI*> vectOfIndicesByIBCType(nbodies);
        vector<E_Int> nPtsPerIBCType(nbodies);//nb de pts projetes sur la surface paroi de type associe 
        for (E_Int notype = 0; notype < nbodies; notype++)
        {
            FldArrayI* indicesByIBCType = new FldArrayI(npts);
            vectOfIndicesByIBCType[notype] = indicesByIBCType;
            nPtsPerIBCType[notype]=0;
        }
        for (E_Int ind = 0; ind < npts; ind++)
        {
            foundw = 0; foundi = 0; 
            E_Int noibctype = -1;
            xc0 = ptrXC[ind]; yc0 = ptrYC[ind]; zc0 = ptrZC[ind];
            /*-----------------------------------------------------------------------------------------------------*/
            /* Step 1: projection  onto the front and bodies following normals -> snearloc resolution for the front*/
            /*-----------------------------------------------------------------------------------------------------*/
            valid = 1;
            dirx0 = ptrNX[ind];
            diry0 = ptrNY[ind];
            dirz0 = ptrNZ[ind];
            //nx = dirx0; ny = diry0; nz = dirz0; //normals towards the outside

#  include "IBC/getIBMPts_projectDirFront.h"

            if (ok == -1) {valid=0; goto projectOrthoBody;}//projectDir on front failed
                        
            //check if projected pts  are close enough - snearloc has been computed more accurately
            delta2 = (xsav-ptrXC[ind])*(xsav-ptrXC[ind])+(ysav-ptrYC[ind])*(ysav-ptrYC[ind])+(zsav-ptrZC[ind])*(zsav-ptrZC[ind]);
            if (delta2 < toldistFactorImage*snearloc)  
            {
                foundi = 1;
                ptrXI[ind] = xsav; ptrYI[ind] = ysav; ptrZI[ind] = zsav;
            }
            else valid=0;

            // projectDir for pt onto the bodies
            dirx0 = sign*dirx0; 
            diry0 = sign*diry0; 
            dirz0 = sign*dirz0;

# include "IBC/getIBMPts_projectDirBodies.h"

            if (ok == -1) {valid = 0; goto projectOrthoBody;}
                        
            delta1 = (xsav-ptrXC[ind])*(xsav-ptrXC[ind])+(ysav-ptrYC[ind])*(ysav-ptrYC[ind])+(zsav-ptrZC[ind])*(zsav-ptrZC[ind]);
            if (delta1 < toldistFactorWall *snearloc)  
            {
                ptrXW[ind] = xsav; ptrYW[ind] = ysav; ptrZW[ind] = zsav;
                foundw = 1;
            }
            else valid = 0;

            if ( valid == 1 ) goto endofpt;//success

            /*--------------------------------------------------------------------------------------*/
            /* Step 2: ortho projection on the bodies + new normals for projectDir on the front     */
            /*--------------------------------------------------------------------------------------*/
            projectOrthoBody:;
            printf("Step1 failed: ortho projection on bodies.\n");

            valid = 1;
            // ortho projection of target pt onto bodies
            pt[0] = xc0; pt[1] = yc0; pt[2] = zc0;
            indp = kdt.getClosest(pt);

# include "IBC/getIBMPts_projectOrthoBodies.h"
            if (ok==-1)  {valid = 0; goto projectOrthoFront;}

            delta1 = (xsav-ptrXC[ind])*(xsav-ptrXC[ind])+(ysav-ptrYC[ind])*(ysav-ptrYC[ind])+(zsav-ptrZC[ind])*(zsav-ptrZC[ind]);
            if ( delta1 < toldistFactorWall *snearloc )  
            {
                ptrXW[ind] = xsav; ptrYW[ind] = ysav; ptrZW[ind] = zsav;
                foundw = 1;
            }
            else valid=0;
        
            // projectDir for pt onto front surfaces using the new direction 
            dirx0 = sign*(xsav-xc0); //sign =1 (internal): follows vector CW
            diry0 = sign*(ysav-yc0); //sign=-1 (external): follows vector WC
            dirz0 = sign*(zsav-zc0);
            nxp = dirx0; nyp = diry0; nzp = dirz0; //normale orientee vers l'exterieur

#  include "IBC/getIBMPts_projectDirFront.h"
            if (ok == -1) { valid=0; goto projectOrthoFront;}
                                 
            delta2 = (xsav-ptrXC[ind])*(xsav-ptrXC[ind])+(ysav-ptrYC[ind])*(ysav-ptrYC[ind])+(zsav-ptrZC[ind])*(zsav-ptrZC[ind]);
            if (delta2 < toldistFactorImage*snearloc)  
            {
                foundi = 1;
                ptrXI[ind] = xsav; ptrYI[ind] = ysav; ptrZI[ind] = zsav;
            }
            else valid = 0;

            if ( valid == 1 ) goto endofpt; // both pts are valid: success
            
            /*---------------------------------------------------------------------------------------------*/
            /* Step 3: ortho projection onto the front + dir onto the obstacle using the ortho direction  */
            /*---------------------------------------------------------------------------------------------*/
            projectOrthoFront:;            
            printf("Step2 failed: ortho projection on front.\n");
            valid = 1;

            // Projection orthogonale du pt cible sur le front
            pt[0] = xc0; pt[1] = yc0; pt[2] = zc0;
            indp = kdtf.getClosest(pt);

# include "IBC/getIBMPts_projectOrthoFront.h"

            if (ok == -1) { valid=0; goto blended; } // nxs non initialise!!
        
            //check if projected pts are close enough
            delta2 = (xsav-ptrXC[ind])*(xsav-ptrXC[ind])+(ysav-ptrYC[ind])*(ysav-ptrYC[ind])+(zsav-ptrZC[ind])*(zsav-ptrZC[ind]);

            if ( delta2 < toldistFactorImage*snearloc )  
            {
                foundi = 1;
                ptrXI[ind] = xsav; ptrYI[ind] = ysav; ptrZI[ind] = zsav;
            }
            else valid = 0;

            //projection to get wall pt
            dirx0 = xc0-xsav; 
            diry0 = yc0-ysav;
            dirz0 = zc0-zsav;
            nxs = -dirx0; nys = -diry0; nzs = -dirz0; //normale orientee vers l'exterieur

# include "IBC/getIBMPts_projectDirBodies.h"

            if ( ok==-1) { valid=0; goto blended;}

            delta1 = (xsav-ptrXC[ind])*(xsav-ptrXC[ind])+(ysav-ptrYC[ind])*(ysav-ptrYC[ind])+(zsav-ptrZC[ind])*(zsav-ptrZC[ind]);
            if ( delta1 < toldistFactorWall *snearloc )  
            {
                ptrXW[ind] = xsav; ptrYW[ind] = ysav; ptrZW[ind] = zsav;
                foundw = 1;
            }
            else valid = 0;

            //check if projected pts  are close enough
            if ( valid == 1 ) goto endofpt; // both pts are valid

            /*------------------------------------------------------------*/
            /*   BLEND OF BOTH DIRECTIONS                                 */
            /*------------------------------------------------------------*/
            blended:;
            printf("Step 3 failed: blending normals...\n");
            valid = 1;

            // projection onto the front
            dirx0 = 0.5*(nxp+nxs);
            diry0 = 0.5*(nyp+nys);
            dirz0 = 0.5*(nzp+nzs);
#  include "IBC/getIBMPts_projectDirFront.h"
            err = 0;
            if (ok == -1) err = 1;
            else 
            {
                delta2 = (xsav-ptrXC[ind])*(xsav-ptrXC[ind])+(ysav-ptrYC[ind])*(ysav-ptrYC[ind])+(zsav-ptrZC[ind])*(zsav-ptrZC[ind]);
                if ( delta2 < toldistFactorImage*snearloc)  
                {
                    foundi = 1;
                    ptrXI[ind] = xsav; ptrYI[ind] = ysav; ptrZI[ind] = zsav;
                }
                else valid = 0;
            }

            //projection onto bodies
            dirx0 = sign*dirx0;
            diry0 = sign*diry0;
            dirz0 = sign*dirz0;

# include "IBC/getIBMPts_projectDirBodies.h"
            if ( ok == -1) err = 1;
            else 
            {
                delta1 = (xsav-ptrXC[ind])*(xsav-ptrXC[ind])+(ysav-ptrYC[ind])*(ysav-ptrYC[ind])+(zsav-ptrZC[ind])*(zsav-ptrZC[ind]);
                if ( delta1 < toldistFactorWall *snearloc )  
                {
                    ptrXW[ind] = xsav; ptrYW[ind] = ysav; ptrZW[ind] = zsav;
                    foundw = 1;
                }
                else valid = 0;

                //check if projected pts  are close enough 
                if (valid == 0) err = 1;
                else goto endofpt;
            }
            if ( err == 1 )
            {
                printf("Warning: all the projections onto the front failed for corrected point %d of zone %d of coordinates (%g,%g,%g).\n", ind, noz, xc0, yc0, zc0);
                pt[0] = xc0; pt[1] = yc0; pt[2] = zc0;
                if ( foundw == 0)
                {
                    printf("Wall point is the closest.\n");
                    indp = kdt.getClosest(pt);
                    ptrXW[ind] = xb2[indp]; ptrYW[ind] = yb2[indp]; ptrZW[ind] = zb2[indp];
                }
                if (foundi == 0)
                {
                    printf("Image point is the closest.\n");
                    indp = kdtf.getClosest(pt);
                    ptrXI[ind] = xf2[indp]; ptrYI[ind] = yf2[indp]; ptrZI[ind] = zf2[indp];
                }
            }
            
            /*------------------------------------------------------------*/
            /*   HOPEFULLY ALL THE CASES ARE VALID NOW */
            /*------------------------------------------------------------*/
            endofpt:;

            E_Int& nptsByType = nPtsPerIBCType[noibctype];
            E_Int* indicesForIBCType = vectOfIndicesByIBCType[noibctype]->begin();
            indicesForIBCType[nptsByType] = ind+1; 
            nptsByType+=1;
        }// ind in zone

        PyObject* PyListIndicesByIBCTypeForZone = PyList_New(0);
        for (E_Int noibctype = 0; noibctype < nbodies; noibctype++)
        {
            E_Int nptsByType = nPtsPerIBCType[noibctype];
            FldArrayI* indicesByIBCTypeL = vectOfIndicesByIBCType[noibctype];
            indicesByIBCTypeL->resize(nptsByType);
            PyObject* tpl0 = K_NUMPY::buildNumpyArray(*indicesByIBCTypeL,1);
            PyList_Append(PyListIndicesByIBCTypeForZone, tpl0); 
            Py_DECREF(tpl0); delete indicesByIBCTypeL;
        }
        PyList_Append(PyListIndicesByIBCType, PyListIndicesByIBCTypeForZone); 
        Py_DECREF(PyListIndicesByIBCTypeForZone);
    }
    delete frontSurfacesPts; delete bodySurfacesPts;

    // Cleaning
    E_Int nboxes = vectOfBodyBoxes.size();
    for (E_Int v0 = 0; v0 < nboxes; v0++)
    {
        vector<BBox3DType*>& boxes = vectOfBodyBoxes[v0];
        E_Int size = boxes.size();
        for (E_Int v = 0; v < size; v++) delete boxes[v];
            delete vectOfBodyBBTrees[v0];
    }
    vectOfBodyBoxes.clear(); vectOfBodyBBTrees.clear();
    nboxes = vectOfFrontBoxes.size();
    for (E_Int v0 = 0; v0 < nboxes; v0++)
    {
        vector<BBox3DType*>& boxes = vectOfFrontBoxes[v0];
        E_Int size = boxes.size();
        for (E_Int v = 0; v < size; v++) delete boxes[v];
            delete vectOfFrontBBTrees[v0];
    }
    vectOfFrontBoxes.clear(); vectOfFrontBBTrees.clear();

    // Sortie
    PyObject* tpl = Py_BuildValue("[OOO]", PyListOfWallPts, PyListOfImagePts, PyListIndicesByIBCType);
    Py_DECREF(PyListOfWallPts);
    Py_DECREF(PyListOfImagePts);
    Py_DECREF(PyListIndicesByIBCType);
    RELEASEZONES; RELEASEBODIES; RELEASEFRONT;
    return tpl;
}
