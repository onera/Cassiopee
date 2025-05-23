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

# include "connector.h"
# include "Nuga/include/BbTree.h"
# include "Nuga/include/KdTree.h"
# include "Nuga/include/ArrayAccessor.h"
# include <math.h>
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
    PyObject *ListOfModelisationHeightsLoc;
    E_Int signOfDist; //if correctedPts are inside bodies: sign = -1, else sign=1 (e.g. for Euler sign=-1, for wall modeling sign=1)
    E_Int depth;//nb of layers of ghost cells
    E_Int isWireModel,isOrthoFirst;
    if (!PYPARSETUPLE_(args, OOOO_ OO_ IIII_,
                       &allCorrectedPts, &ListOfSnearsLoc, &ListOfModelisationHeightsLoc, &bodySurfaces, &frontSurfaces, 
                       &normalNames, &signOfDist, &depth, &isWireModel, &isOrthoFirst)) return NULL;

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

    // extract list of ModelisationHeightsLoc
    if (PyList_Size(ListOfModelisationHeightsLoc) == 0)
    {
        PyErr_SetString(PyExc_TypeError, 
                        "getIBMPtsWithFront: 3rd argument is an empty list.");
        return NULL;
    }

    vector<E_Float> vectOfModelisationHeightsLoc;

    for (int i = 0; i < nsnears; i++)
    {
        tpl0 = PyList_GetItem(ListOfModelisationHeightsLoc, i);
        if (PyFloat_Check(tpl0) == 0)
        {
            PyErr_SetString(PyExc_TypeError, 
                            "getIBMPtsWithFront: not a valid value for modelisation height.");
            return NULL;
        } 
        else vectOfModelisationHeightsLoc.push_back(PyFloat_AsDouble(tpl0));
    }
  
    E_Int sign = -signOfDist; // sens de projection sur la paroi
    if (isWireModel==1){
      sign = signOfDist;
    }
    //tolerance of distance between corrected pts and wall/image to control projection
    //distance max for image pt to its corrected pt : (depth+1)*sqrt(2)*snearloc 
    E_Float toldistFactorImage;
    if (signOfDist==-1) // interior points:  we move away the front from one additional layer
        toldistFactorImage = ((depth+1)*1.1)*((depth+1)*1.1)*3.;
    else 
        toldistFactorImage = (depth*1.1)*(depth*1.1)*3.;
    //distance max for wall pt to its corrected pt ( depth)*sqrt(3)*snearloc
    E_Float toldistFactorWall = (depth*1.1)*(depth*1.1)*3.;

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
        PyObject* l = PyList_GetItem(normalNames, v);
        if (PyString_Check(l))
        {
            var = PyString_AsString(l);
            varsn.push_back(var);
        }
#if PY_VERSION_HEX >= 0x03000000
        else if (PyUnicode_Check(l)) 
        {
            var = (char*)PyUnicode_AsUTF8(l);
            varsn.push_back(var);
        } 
#endif
        else
        {
            PyErr_SetString(PyExc_TypeError,
                            "getIBMPts: invalid string for normal component.");
            return NULL;
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
    PyObject* PyListOfProjectionType = PyList_New(0);
    
    E_Int nvarOut = 3;
    char varStringOut[K_ARRAY::VARSTRINGLENGTH]; 
    strcpy(varStringOut,"CoordinateX,CoordinateY,CoordinateZ");
    char eltTypeOut[8]; strcpy(eltTypeOut,"NODE");
    FldArrayI cnl(0);
    E_Int nelts = 0;
    vector<E_Float*> xit; vector<E_Float*> yit; vector<E_Float*> zit;
    vector<E_Float*> xwt; vector<E_Float*> ywt; vector<E_Float*> zwt;
    vector<E_Float*> xPTt; 

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

	PyObject* tplPT = K_ARRAY::buildArray(1, "ProjectionType", npts, nelts, -1, eltTypeOut, false);
        E_Float* coordPTp = K_ARRAY::getFieldPtr(tplPT);
        FldArrayF coordPT(npts, 1, coordPTp, true); // non initialized yet
        coordPT.setOneField(*unstrF[no],posxt[no],1);
        xPTt.push_back(coordPT.begin(1));
        PyList_Append(PyListOfProjectionType, tplPT); Py_DECREF(tplPT);

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
    E_Float dirx0, diry0, dirz0, xsf, ysf, zsf, xsb, ysb, zsb;
    E_Float xc0, yc0, zc0, xw0, yw0, zw0, xi0, yi0, zi0;
    E_Float dist2, distl;
    vector<E_Int> indicesBB; 
    E_Float pr1[3]; E_Float pr2[3]; E_Float pt[3];
    E_Float p0[3]; E_Float p1[3]; E_Float p2[3]; E_Float p[3];
    E_Int oriented = 1;
    E_Int ok, notri, indvert1, indvert2, indvert3, indp;
    E_Float rx, ry, rz, rad; 

    //distance of corrected pts to wall pts and image pts 
    //E_Float delta1, delta2;
    //Corresponding directions
    E_Float nxp, nyp, nzp, nxs, nys, nzs;
    E_Int *cnVert1, *cnVert2, *cnVert3;
    E_Float edgeLen1, edgeLen2;
    PyObject* PyListIndicesByIBCType = PyList_New(0);
    E_Float xf_ortho, yf_ortho, zf_ortho, xb_ortho, yb_ortho, zb_ortho;
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
        E_Float* ptrPT = xPTt[noz];
        E_Float snearloc = vectOfSnearsLoc[noz];
        E_Float heightloc = vectOfModelisationHeightsLoc[noz];
        //distance max for image pt to its corrected pt : height*sqrt(3)

        //old way of setting max tolerances
        heightloc = heightloc*1.1 + 3*snearloc*sqrt(3.); // for 2nd image point
        heightloc = heightloc*heightloc;

        snearloc = snearloc + 3*snearloc*sqrt(3); // for 2nd image point
        snearloc = snearloc*snearloc;

        E_Float distMaxF2 = max(toldistFactorImage*snearloc, heightloc);// distance au carre maximale des pts cibles au front via depth ou modelisationHeight
        E_Float distMaxB2 = max(toldistFactorWall*snearloc, heightloc);// distance au carre maximale des pts cibles au projete paroi via depth ou modelisationHeight

        //new way of setting max tolerances
        // heightloc = heightloc*heightloc;
        // snearloc = snearloc*snearloc;
        // E_Float distMaxF2 = max(toldistFactorImage*snearloc, 4*heightloc); // squared maximum projection distance for target points based on local near-wall resolution or modeling height
        // E_Float distMaxB2 = max(toldistFactorWall*snearloc, 4*heightloc); // squared maximum projection distance for target points based on local near-wall resolution or modeling height
        // These max distances are based on 2*hmod instead of hmod to allow some tolerance for complex geometries

        vector<FldArrayI*> vectOfIndicesByIBCType(nbodies);
        vector<E_Int> nPtsPerIBCType(nbodies);//nb de pts projetes sur la surface paroi de type associe 
        for (E_Int notype = 0; notype < nbodies; notype++)
        {
            FldArrayI* indicesByIBCType = new FldArrayI(npts);
            vectOfIndicesByIBCType[notype] = indicesByIBCType;
            nPtsPerIBCType[notype] = 0;
        }
        FldArrayI typeProj(npts); typeProj.setAllValuesAtNull();
        E_Int nType1=0, nType2=0, nType3=0, nType4=0;
        for (E_Int ind = 0; ind < npts; ind++)
        {
            E_Int noibctype = -1;
            xc0 = ptrXC[ind]; yc0 = ptrYC[ind]; zc0 = ptrZC[ind];
            E_Float distF1 = -1.; E_Float distB1 = -1.;
            E_Float distF2 = -1.; E_Float distB2 = -1.;
            E_Float distF3 = -1.; E_Float distB3 = -1.;
	    
            /*-----------------------------------------------------------------------------------------------------*/
            /* Step 1: projection  onto the front and bodies following normals -> snearloc resolution for the front*/
            /*-----------------------------------------------------------------------------------------------------*/
            E_Int found = 0;

            dirx0 = ptrNX[ind]; diry0 = ptrNY[ind]; dirz0 = ptrNZ[ind];
            if (isOrthoFirst == 1) {goto ortho_projection;}
#  include "IBC/getIBMPts_projectDirFront.h"

            if (ok > -1) // projection found
            {
                distF1 = (xsf-xc0)*(xsf-xc0)+(ysf-yc0)*(ysf-yc0)+(zsf-zc0)*(zsf-zc0);            

                // projectDir for pt onto the bodies
                dirx0 = sign*dirx0; 
                diry0 = sign*diry0; 
                dirz0 = sign*dirz0;

                # include "IBC/getIBMPts_projectDirBodies.h"
                if (ok > -1) 
                {
                    distB1 = (xsb-xc0)*(xsb-xc0)+(ysb-yc0)*(ysb-yc0)+(zsb-zc0)*(zsb-zc0);

                    //check distance
                    if (distB1 <= distMaxB2 && distF1 <= distMaxF2)
                    {
                        found = 1;
                        ptrXW[ind] = xsb; ptrYW[ind] = ysb; ptrZW[ind] = zsb;
                        ptrXI[ind] = xsf; ptrYI[ind] = ysf; ptrZI[ind] = zsf;
                        typeProj[ind]=1; nType1+=1;
                        goto end;
                    }
                }
            }

            // found = 1: image and wall IBM points have been found and set
	ortho_projection:
            if (found == 0)
            {
                /*--------------------------------------------------------------------------------------*/
                /* Step 2: ortho projection on the bodies + new normals for projectDir on the front     */
                /*--------------------------------------------------------------------------------------*/
                // ortho projection of target pt onto bodies
                pt[0] = xc0; pt[1] = yc0; pt[2] = zc0;
                indp = kdt.getClosest(pt);
                // std::cout << "found = 0" << std::endl;
                # include "IBC/getIBMPts_projectOrthoBodies.h"
		
                if (ok == 1)
                { 
                    distB2 = (xsb-xc0)*(xsb-xc0)+(ysb-yc0)*(ysb-yc0)+(zsb-zc0)*(zsb-zc0);

                    xb_ortho = xsb; yb_ortho = ysb; zb_ortho=zsb;

                    // projectDir for pt onto front surfaces using the new direction
                    dirx0 = sign*(xsb-xc0); //sign = 1 (internal): follows vector CW
                    diry0 = sign*(ysb-yc0); //sign =-1 (external): follows vector WC
                    dirz0 = sign*(zsb-zc0);
	    
                    nxp = dirx0; nyp = diry0; nzp = dirz0; //normale orientee vers l'exterieur
                    #  include "IBC/getIBMPts_projectDirFront.h"
		    if (ok == 1)
                    {
                        distF2 = (xsf-xc0)*(xsf-xc0)+(ysf-yc0)*(ysf-yc0)+(zsf-zc0)*(zsf-zc0);
                        if (distF2  <= distMaxF2 && distB2 <= distMaxB2)
                        {
                            ptrXW[ind] = xsb; ptrYW[ind] = ysb; ptrZW[ind] = zsb;
                            ptrXI[ind] = xsf; ptrYI[ind] = ysf; ptrZI[ind] = zsf;
                            typeProj[ind]=2; nType2+=1;
			    goto end;
                        }
                    }
                }
                else {printf("WARNING: getIBMPtsWithFront: projectOrthoBodies FAILED !!!\n");}
                /*---------------------------------------------------------------------------------------------*/
                /* Step 3: ortho projection onto the front + dir onto the obstacle using the ortho direction  */
                /*---------------------------------------------------------------------------------------------*/
                // std::cout << "step 2 failed" << std::endl;
                pt[0] = xc0; pt[1] = yc0; pt[2] = zc0;
                indp = kdtf.getClosest(pt);
                # include "IBC/getIBMPts_projectOrthoFront.h"
                if (ok == 1)
                {
                    distF3 = (xsf-xc0)*(xsf-xc0)+(ysf-yc0)*(ysf-yc0)+(zsf-zc0)*(zsf-zc0);
                    xf_ortho=xsf; yf_ortho=ysf; zf_ortho=zsf;

                    //projectDir to get wall pt
                    dirx0 = xc0-xsf; diry0 = yc0-ysf; dirz0 = zc0-zsf;
                    nxs = -dirx0; nys = -diry0; nzs = -dirz0; //normale orientee vers l'exterieur
                    # include "IBC/getIBMPts_projectDirBodies.h" 
                    if (ok == 1)
                    {
                        distB3 = (xsb-xc0)*(xsb-xc0)+(ysb-yc0)*(ysb-yc0)+(zsb-zc0)*(zsb-zc0);
                        if (distF3 <= distMaxF2 && distB3 <= distMaxB2) 
                        {
                            ptrXW[ind] = xsb; ptrYW[ind] = ysb; ptrZW[ind] = zsb;
                            ptrXI[ind] = xsf; ptrYI[ind] = ysf; ptrZI[ind] = zsf;
                            typeProj[ind]=3; nType3+=1;
                            goto end;
                        }
                    }
                }
                else {printf("WARNING: getIBMPtsWithFront: projectOrthoFront FAILED !!!\n");}

               /*---------------------------------------------------------------------------------------------*/
               /* Step  4:  no valid projection found. Minimum distance chosen (projectOrtho)*/
               /*---------------------------------------------------------------------------------------------*/
                // std::cout << "step 3 failed" << std::endl;
                // printf("Warning: minimum distance chosen for point of indice %d of zone %d.\n",ind,noz);
                // printf(" xc = %g,%g,%g \n", xc0, yc0, zc0);
                // printf(" xw = %g,%g,%g \n", xb_ortho, yb_ortho, zb_ortho);
                // printf(" xi = %g,%g,%g \n", xf_ortho, yf_ortho, zf_ortho);
                typeProj[ind]=4; nType4+=1;

                ptrXW[ind] = xb_ortho; ptrYW[ind] = yb_ortho; ptrZW[ind] = zb_ortho;
                ptrXI[ind] = xf_ortho; ptrYI[ind] = yf_ortho; ptrZI[ind] = zf_ortho;
            }// found CAS 1 = 0         
            end:;

	    ptrPT[ind]=typeProj[ind];
	    
            E_Int& nptsByType = nPtsPerIBCType[noibctype];
            E_Int* indicesForIBCType = vectOfIndicesByIBCType[noibctype]->begin();
            indicesForIBCType[nptsByType] = ind+1; 
            nptsByType += 1;
            // BILAN 
        }// ind in zone
        // printf(" ZONE %d Nb de type 1 : %d, type2 = %d, type3 = %d, type4 = %d\n", noz, nType1, nType2,nType3,nType4);
        if (nType3 > 0 || nType4 > 0)
        {
            printf("Warning getIBMPtsWithFront: #Zone " SF_D_ " has " SF_D_ " pts of type 3 or 4 ! (Bilan: t1/t2/t3/t4 = %.2f%%/%.2f%%/%.2f%%/%.2f%%)\n", noz, nType2+nType3+nType4,
            nType1/E_Float(nType1+nType2+nType3+nType4)*100.,
            nType2/E_Float(nType1+nType2+nType3+nType4)*100.,
            nType3/E_Float(nType1+nType2+nType3+nType4)*100.,
            nType4/E_Float(nType1+nType2+nType3+nType4)*100. );
        }

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
    PyObject* tpl = Py_BuildValue("[OOOO]", PyListOfWallPts, PyListOfImagePts, PyListIndicesByIBCType, PyListOfProjectionType);
    Py_DECREF(PyListOfWallPts);
    Py_DECREF(PyListOfImagePts);
    Py_DECREF(PyListOfProjectionType);
    Py_DECREF(PyListIndicesByIBCType);
    RELEASEZONES; RELEASEBODIES; RELEASEFRONT;
    return tpl;
}
