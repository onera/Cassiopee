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

// ============================================================================
/* Get the wall and image points if front is not provided 
    wall pts are obtained by projDir (projOrtho, closest pt if impossible)
    image pts are obtained by translation of corrected pts of hi/he*/
// ============================================================================
PyObject* K_CONNECTOR::getIBMPtsWithoutFront(PyObject* self, PyObject* args)
{
    PyObject *allCorrectedPts, *bodySurfaces, *normalNames, *distName;
    E_Int signOfDist; //if correctedPts are inside bodies: sign = -1, else sign=1
    if (!PYPARSETUPLE_(args, OOOO_ I_,
                       &allCorrectedPts, &bodySurfaces, 
                       &normalNames, &distName, &signOfDist))
        return NULL;

    E_Int sign = -signOfDist; // sens de projection sur la paroi

    // check distname
    char* distname;
    if (PyString_Check(distName)) distname = PyString_AsString(distName);
#if PY_VERSION_HEX >= 0x03000000
    else if (PyUnicode_Check(distName)) distname = (char*)PyUnicode_AsUTF8(distName);
#endif
    else
    {    
        PyErr_SetString(PyExc_TypeError, 
                        "getIBMPtsWithoutFront: distName must be a string.");
        return NULL;
    }

    // Check normal components
    if (PyList_Check(normalNames) == 0)
    {
        PyErr_SetString(PyExc_TypeError, 
                        "getIBMPtsWithoutFront: normal vars must be a list of strings.");
        return NULL;
    }
    E_Int nvars = PyList_Size(normalNames);
    if (nvars != 3)
    {
        PyErr_SetString(PyExc_TypeError, 
                        "getIBMPtsWithoutFront: normal vars must be a 3-component vector.");
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
                            "getIBMPtsWithoutFront: invalid string for normal component.");
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
    E_Bool skipNoCoord = true;
    E_Bool skipStructured = true;
    E_Bool skipUnstructured = false;
    E_Bool skipDiffVars = false;
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
    vector<E_Int> posnxt; vector<E_Int> posnyt; vector<E_Int> posnzt; vector<E_Int> posdist;
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
                
        posx1++; posy1++; posz1++;
        posnxt.push_back(posx1); posnyt.push_back(posy1); posnzt.push_back(posz1);
        
        posx1 = K_ARRAY::isNamePresent(distname, unstrVarString[no]);   
        if ( posx1 == -1)
        {
            PyErr_SetString(PyExc_TypeError,"getIBMPts: dist variable is not present in 1st arg.");
            RELEASEZONES; return NULL;
        }
        posx1++; posdist.push_back(posx1);
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
        PyErr_SetString(PyExc_TypeError,"getIBMPts: 2nd arg is not valid.");
        RELEASEZONES;
        RELEASEBODIES;
        return NULL;
    }
    for (E_Int no = 0; no < nbodies; no++)
    {
        if (K_STRING::cmp(eltTypeb[no], "TRI") != 0) 
        {
            PyErr_SetString(PyExc_TypeError,"getIBMPts: 2nd arg must be a TRI type zone.");
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
    for (E_Int v = 0; v < nbodies; v++)
        nptsBodiesMax += unstrbF[v]->getSize();
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

    /*--------------------------------------------------------*/
    // projectDir of all the points onto the bodies-> wall pts
    /*--------------------------------------------------------*/    
    E_Float tol = K_CONST::E_GEOM_CUTOFF;
    E_Float dirn, dirx0, diry0, dirz0;
    E_Float xsb, ysb, zsb;
    E_Float dist0, xc0, yc0, zc0, xw0, yw0, zw0, xi0, yi0, zi0;
    E_Float dist2, distl;
    vector<E_Int> indicesBB; 
    E_Float pr1[3]; E_Float pr2[3]; E_Float pt[3];
    E_Float p0[3]; E_Float p1[3]; E_Float p2[3]; E_Float p[3];
    E_Int oriented = 1;
    E_Int ok, notri, indp;
    E_Float rx, ry, rz, rad;
    // Corresponding directions
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
        E_Float* distp  = correctedPts->begin(posdist[noz]);

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
            xc0 = ptrXC[ind]; yc0 = ptrYC[ind]; zc0 = ptrZC[ind];
            E_Int noibctype = -1;

            /*--------------------------------------------------*/
            /*  STEP 1/ projection following normals onto bodies*/
            /*--------------------------------------------------*/
            dirx0 = sign*ptrNX[ind]; 
            diry0 = sign*ptrNY[ind]; 
            dirz0 = sign*ptrNZ[ind];

# include "IBC/getIBMPts_projectDirBodies.h"

            if ( ok == -1) goto projectOrthoBody;
            else 
            {
                ptrXW[ind] = xsb; ptrYW[ind] = ysb; ptrZW[ind] = zsb;
                goto endofwpt;
            }
            /*----------------------------------------*/
            /* STEP2/ ortho projection on the bodies  */
            /*----------------------------------------*/
            projectOrthoBody:;            

            // Projection orthogonale du pt cible sur la paroi
            pt[0] = xc0; pt[1] = yc0; pt[2] = zc0;
            indp = kdt.getClosest(pt);

# include "IBC/getIBMPts_projectOrthoBodies.h"

            if (ok==-1) //closest pt
            {
                xsb = xb2[indp]; ysb = yb2[indp]; zsb = zb2[indp];
            }

            ptrXW[ind] = xsb; ptrYW[ind] = ysb; ptrZW[ind] = zsb;
            endofwpt:;            
            /*----------------------------------------*/
            /* STEP3/ determination of image pts      */
            /*----------------------------------------*/
            dirx0 = (ptrXC[ind]-ptrXW[ind]);
            diry0 = (ptrYC[ind]-ptrYW[ind]);
            dirz0 = (ptrZC[ind]-ptrZW[ind]);
            dirn = sqrt(dirx0*dirx0+diry0*diry0+dirz0*dirz0);
            dist0 = distp[ind]* signOfDist/dirn;
            ptrXI[ind] = ptrXW[ind] + dirx0*dist0;
            ptrYI[ind] = ptrYW[ind] + diry0*dist0;
            ptrZI[ind] = ptrZW[ind] + dirz0*dist0;

            E_Int& nptsByType = nPtsPerIBCType[noibctype];
            E_Int* indicesForIBCType = vectOfIndicesByIBCType[noibctype]->begin();
            indicesForIBCType[nptsByType] = ind+1; 
            nptsByType+=1;
        }//ind

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
    delete bodySurfacesPts;
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


    // Sortie
    PyObject* tpl = Py_BuildValue("[OOO]", PyListOfWallPts, PyListOfImagePts, PyListIndicesByIBCType);
    Py_DECREF(PyListOfWallPts);
    Py_DECREF(PyListOfImagePts);
    Py_DECREF(PyListIndicesByIBCType);
    RELEASEZONES; RELEASEBODIES;
    return tpl;
}
