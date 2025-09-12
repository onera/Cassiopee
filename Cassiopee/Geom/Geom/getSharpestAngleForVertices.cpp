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

# include "geom.h"
using namespace K_FLD;
using namespace std;

//=============================================================================
/* Compute the sharpest angle alpha between elements of a surface. 
   The angle is returned at vertices */
//=============================================================================
PyObject* K_GEOM::getSharpestAngleForVertices(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float dirVect[3];
  dirVect[0] = 0.;  dirVect[1] = 0.;  dirVect[2] = 1.;
  if (!PYPARSETUPLE_(args, O_, &array)) return NULL; 
  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, 
                                     ni, nj, nk, cn, eltType);
 
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError, "getSharpestAngleForVertices: coordinates not found in array.");
    RELEASESHAREDU(array,f, cn); return NULL;
  }
  posx++; posy++; posz++;

  if (res == 1) 
  { 
    if (nj > 1 || nk > 1) 
    {
      PyErr_SetString(PyExc_TypeError, "getSharpestAngleForVertices: structured array must be 1D.");
      RELEASESHAREDS(array,f); return NULL;
    }
    E_Int npts = f->getSize();
    PyObject* tpl = K_ARRAY::buildArray(1, "alpha", npts, 1, 1);
    E_Float* alpt = K_ARRAY::getFieldPtr(tpl);
    E_Float ptA1[3]; E_Float ptB1[3];  E_Float ptA2[3]; E_Float ptB2[3]; 
    E_Float* xt = f->begin(posx); E_Float* yt = f->begin(posy); E_Float* zt = f->begin(posz);
    E_Int i, im, ip;
    E_Float alpha0;
    alpt[0] = 0.; alpt[npts-1] = 0.; 
    E_Int closed = 0;
    if (K_FUNC::fEqualZero(xt[0]-xt[npts-1]) == true && 
        K_FUNC::fEqualZero(yt[0]-yt[npts-1]) == true && 
        K_FUNC::fEqualZero(zt[0]-zt[npts-1]) == true) closed = 1;
    for (i = 1; i < npts-1; i++)
    {
      im = i-1; ip = i+1;
      ptA1[0] = xt[i]; ptA1[1] = yt[i]; ptA1[2] = zt[i];
      ptB1[0] = xt[im]; ptB1[1] = yt[im]; ptB1[2] = zt[im];
      ptA2[0] = xt[i]; ptA2[1] = yt[i]; ptA2[2] = zt[i];
      ptB2[0] = xt[ip]; ptB2[1] = yt[ip]; ptB2[2] = zt[ip];
      alpha0 = K_COMPGEOM::getAlphaAngleBetweenBars(ptA1, ptB1, ptA2, ptB2, dirVect);
      if (alpha0 != -1000.) alpt[i] = alpha0;
      else alpt[i] = 0.;
    }
    if (closed == 1) // loop
    {
      i = 0; ip = 1; im = npts-2;
      ptA1[0] = xt[i]; ptA1[1] = yt[i]; ptA1[2] = zt[i];
      ptB1[0] = xt[im]; ptB1[1] = yt[im]; ptB1[2] = zt[im];
      ptA2[0] = xt[i]; ptA2[1] = yt[i]; ptA2[2] = zt[i];
      ptB2[0] = xt[ip]; ptB2[1] = yt[ip]; ptB2[2] = zt[ip];
      alpha0 = K_COMPGEOM::getAlphaAngleBetweenBars(ptA1, ptB1, ptA2, ptB2, dirVect);
      if (alpha0 != -1000.) {alpt[0] = alpha0; alpt[npts-1] = alpha0;}
    }
    RELEASESHAREDS(array,f); 
    return tpl;
  }
  else // ( res == 2 ) 
  {
    if (strcmp(eltType, "NGON") == 0) // NGON
    {
      E_Int npts = f->getSize();
      E_Int* cnp = cn->begin(); // pointeur sur la connectivite NGon
      E_Int sizeFN = cnp[1]; //  taille de la connectivite Face/Noeuds
      E_Int nelts = cnp[sizeFN+2];  // nombre total d elements
      PyObject* tpl = K_ARRAY::buildArray(1, "alpha", npts, nelts, -1, eltType, 
                                          false, cn->getSize());
      E_Float* alpt = K_ARRAY::getFieldPtr(tpl);
      E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
      FldArrayI cnn(cn->getSize(), 1, cnnp, true); cnn = *cn;
      E_Int nfaces = cnp[0];
      FldArrayI pos; K_CONNECT::getPosElts(*cn, pos);
      FldArrayI posFaces(nfaces); K_CONNECT::getPosFaces(*cn, posFaces);
      //E_Int* posFacesp = posFaces.begin();
      vector< vector<E_Int> > cVF(npts); K_CONNECT::connectNG2VF(*cn, cVF);
      FldArrayI cFE; K_CONNECT::connectNG2FE(*cn, cFE);
      FldArrayI dimElts(npts); K_CONNECT::getDimElts(*cn, dimElts);
      E_Int dim = dimElts[0]; // dimension de la zone
      E_Int indf, e1, e2, ind1, ind2, indl;
      E_Int* pt; E_Float dalphamin, alpha;
      E_Float xbf, ybf, zbf, alpha0;
      E_Float ptA[3]; E_Float ptB[3]; E_Float ptC[3]; E_Float ptD[3];
      E_Float* xt = f->begin(posx); E_Float* yt = f->begin(posy); E_Float* zt = f->begin(posz);
      vector<E_Int> indices; E_Int c;

      if (dim == 1)
      {
        PyErr_SetString(PyExc_TypeError, "getSharpestAngleForVertices: array must be NGON2D.");
        RELEASESHAREDU(array,f, cn); return NULL;       
      }
      else if (dim == 2)
      {
        for (E_Int ind = 0; ind < npts; ind++)
        {
          dalphamin = 360.; c =0; alpha = 0.;

          // Get the faces connected to ind
          vector<E_Int>& f = cVF[ind];
          E_Int nf = f.size();
          for (E_Int k = 0; k < nf; k++)
          {
            indf = f[k]; // indice face
            pt = cnp+posFaces[indf-1];
            ind1 = pt[1]-1; ind2 = pt[2]-1;
            ptA[0] = xt[ind1]; ptA[1] = yt[ind1]; ptA[2] = zt[ind1];
            ptB[0] = xt[ind2]; ptB[1] = yt[ind2]; ptB[2] = zt[ind2];

            e1 = cFE(indf-1, 1); // element gauche de la face
            e2 = cFE(indf-1, 2);
            if (e1 > 0 && e2 > 0)
            {
              K_CONNECT::getVertexIndices(cn->begin(), posFaces.begin(), pos[e1-1], indices);
              E_Int nv = indices.size();
              xbf = 0.; ybf = 0.; zbf = 0.;
              for (E_Int m = 0; m < nv; m++)
              {
                indl = indices[m]-1;
                xbf += xt[indl]; ybf += yt[indl]; zbf += zt[indl];
              }
              ptC[0] = xbf/nv; ptC[1] = ybf/nv; ptC[2] = zbf/nv;
              
              K_CONNECT::getVertexIndices(cn->begin(), posFaces.begin(), pos[e2-1], indices);
              nv = indices.size();
              xbf = 0.; ybf = 0.; zbf = 0.;
              for (E_Int m = 0; m < nv; m++)
              {
                indl = indices[m]-1;
                xbf += xt[indl]; ybf += yt[indl]; zbf += zt[indl]; 
              }
              ptD[0] = xbf/nv; ptD[1] = ybf/nv; ptD[2] = zbf/nv;
              alpha0 = K_COMPGEOM::getAlphaAngleBetweenTriangles(ptA, ptB, ptC, ptA, ptD, ptB);
              if (alpha0 != -1000.)
              {
                c++;
                if (alpha0 < 180. && K_FUNC::E_abs(alpha0-90.) < dalphamin) {dalphamin = K_FUNC::E_abs(alpha0-90.); alpha = alpha0;}
                else if (alpha0 >= 180. && K_FUNC::E_abs(alpha0-270.) < dalphamin) {dalphamin = K_FUNC::E_abs(alpha0-270.); alpha = alpha0;}
              }
            }
          }
          if (c > 0) alpt[ind] = alpha;
          else alpt[ind] = 0.;
        }
        RELEASESHAREDU(array,f, cn); return tpl;
      }
      else // dim = 3
      {
        // delete tpl
        PyErr_SetString(PyExc_TypeError, "getSharpestAngleForVertices: array must be NGON2D.");
        RELEASESHAREDU(array,f, cn); return NULL;        
      }
      
    }
    else // BASIC ELTS
    {
      if (strcmp(eltType, "TRI")  != 0 && 
          strcmp(eltType, "QUAD") != 0 &&
          strcmp(eltType, "BAR")  != 0)
      {
        PyErr_SetString(PyExc_TypeError, "getSharpestAngleForVertices: array must be TRI or QUAD.");
        RELEASESHAREDU(array,f, cn); return NULL;
      }

      E_Int npts = f->getSize();
      E_Int type = cn->getNfld(); E_Int nelts = cn->getSize();
      PyObject* tpl = K_ARRAY::buildArray(1, "alpha", npts, nelts, -1, eltType);
      E_Float* alpt = K_ARRAY::getFieldPtr(tpl);
      E_Int* cnp = K_ARRAY::getConnectPtr(tpl);
      FldArrayI cno(nelts, type, cnp, true); cno = *cn;

      E_Float ptA1[3]; E_Float ptB1[3]; E_Float ptC1[3]; E_Float ptD1[3];
      E_Float ptA2[3]; E_Float ptB2[3]; E_Float ptC2[3]; E_Float ptD2[3];
      vector< vector<E_Int> > cVE(npts);
      K_CONNECT::connectEV2VE(cno, cVE);
      E_Float* xt = f->begin(posx); E_Float* yt = f->begin(posy); E_Float* zt = f->begin(posz);
      E_Int indA1, indB1, indC1, indD1, indA2, indB2, indC2, indD2;
      E_Int et1, et2, c; E_Float alpha0, alpha;
      E_Int* cn1 = cno.begin(1); E_Int* cn2 = cno.begin(2);
      E_Float dalphamin = 360.;
      if (type == 3) // TRI
      {
        E_Int* cn3 = cno.begin(3);
        for (E_Int ind = 0; ind < npts; ind++)
        {
          vector<E_Int>& elts = cVE[ind];
          E_Int neltsv = elts.size();
          c = 0; alpha = 0.; dalphamin = 360.;
          for (E_Int noet1 = 0; noet1 < neltsv; noet1++)
          {
            et1 = elts[noet1];
            indA1 = cn1[et1]-1; indB1 = cn2[et1]-1; indC1 = cn3[et1]-1;
            ptA1[0] = xt[indA1]; ptA1[1] = yt[indA1]; ptA1[2] = zt[indA1];
            ptB1[0] = xt[indB1]; ptB1[1] = yt[indB1]; ptB1[2] = zt[indB1];
            ptC1[0] = xt[indC1]; ptC1[1] = yt[indC1]; ptC1[2] = zt[indC1];
            for (E_Int noet2 = noet1+1; noet2 < neltsv; noet2++)
            {
              et2 = elts[noet2];
              indA2 = cn1[et2]-1; indB2 = cn2[et2]-1; indC2 = cn3[et2]-1;
              ptA2[0] = xt[indA2]; ptA2[1] = yt[indA2]; ptA2[2] = zt[indA2];
              ptB2[0] = xt[indB2]; ptB2[1] = yt[indB2]; ptB2[2] = zt[indB2];
              ptC2[0] = xt[indC2]; ptC2[1] = yt[indC2]; ptC2[2] = zt[indC2];                  
              alpha0 = K_COMPGEOM::getAlphaAngleBetweenTriangles(ptA1, ptB1, ptC1, ptA2, ptB2, ptC2);
              if (alpha0 != -1000.)
              {
                c++;
                if (alpha0 < 180. && K_FUNC::E_abs(alpha0-90.) < dalphamin) {dalphamin = K_FUNC::E_abs(alpha0-90.); alpha = alpha0;}
                else if (alpha0 >= 180. && K_FUNC::E_abs(alpha0-270.) < dalphamin) {dalphamin = K_FUNC::E_abs(alpha0-270.); alpha = alpha0;}
              }
            }
          }
          if ( c > 0 ) alpt[ind] = alpha;
          else alpt[ind] = 0.;
        }
        RELEASESHAREDU(array,f, cn); return tpl;
      }
      else if (type == 4)  //QUAD
      {
        E_Int* cn3 = cno.begin(3); E_Int* cn4 = cno.begin(4);
        for (E_Int ind = 0; ind < npts; ind++)
        {
          vector<E_Int>& elts = cVE[ind];
          E_Int neltsv = elts.size();
          c = 0; alpha = 0.; dalphamin = 360.;
          for (E_Int noet1 = 0; noet1 < neltsv; noet1++)
          {
            et1 = elts[noet1];
            indA1 = cn1[et1]-1; indB1 = cn2[et1]-1; indC1 = cn3[et1]-1; indD1 = cn4[et1]-1;
            ptA1[0] = xt[indA1]; ptA1[1] = yt[indA1]; ptA1[2] = zt[indA1];
            ptB1[0] = xt[indB1]; ptB1[1] = yt[indB1]; ptB1[2] = zt[indB1];
            ptC1[0] = xt[indC1]; ptC1[1] = yt[indC1]; ptC1[2] = zt[indC1];
            ptD1[0] = xt[indD1]; ptD1[1] = yt[indD1]; ptD1[2] = zt[indD1];
            for (E_Int noet2 = noet1+1; noet2 < neltsv; noet2++)
            {
              et2 = elts[noet2];
              indA2 = cn1[et2]-1; indB2 = cn2[et2]-1; indC2 = cn3[et2]-1; indD2 = cn4[et2]-1;
              ptA2[0] = xt[indA2]; ptA2[1] = yt[indA2]; ptA2[2] = zt[indA2];
              ptB2[0] = xt[indB2]; ptB2[1] = yt[indB2]; ptB2[2] = zt[indB2];
              ptC2[0] = xt[indC2]; ptC2[1] = yt[indC2]; ptC2[2] = zt[indC2];  
              ptD2[0] = xt[indD2]; ptD2[1] = yt[indD2]; ptD2[2] = zt[indD2];                
              alpha0 = K_COMPGEOM::getAlphaAngleBetweenQuads(ptA1, ptB1, ptC1, ptD1, 
                                                             ptA2, ptB2, ptC2, ptD2);
              if (alpha0 != -1000.)
              {
                c++;
                if ( alpha0 < 180. && K_FUNC::E_abs(alpha0-90.) < dalphamin ) {dalphamin = K_FUNC::E_abs(alpha0-90.); alpha = alpha0;}
                else if ( alpha0 >= 180. && K_FUNC::E_abs(alpha0-270.) < dalphamin ) {dalphamin = K_FUNC::E_abs(alpha0-270.); alpha = alpha0;}
              }
            }
          }
          if ( c > 0 ) { alpt[ind] = alpha; }
          else alpt[ind] = 0.;
        }

        RELEASESHAREDU(array,f, cn); return tpl;  
      }
      else //BAR
      {
        E_Int indm, indp, e1, e2;
        for (E_Int ind = 0; ind < npts; ind++)
        {       
          vector<E_Int>& elts = cVE[ind];
          E_Int neltsv = elts.size();
          if ( neltsv == 2 ) 
          {
            e1 = elts[0];
            e2 = elts[1];
            indA1 = cn1[e1]-1;
            indB1 = cn2[e1]-1;
            indA2 = cn1[e2]-1;
            indB2 = cn2[e2]-1;
            if (ind == indA1 && ind == indB2)
            {
              indp = indB1;
              indm = indA2;
            }
            else if (ind == indB1 && ind == indA2)
            {
              indp = indB2;
              indm = indA1;
            }
            else if (ind == indA1 && ind == indA2)
            {
              // reverse case : ambiguous
              indp = indB1;
              indm = indB2;
            }
            else
            {
              // reverse case : ambiguous
              indp = indA1;
              indm = indA2;
            }           
          }
          else if ( neltsv == 1)
          {
            e1 = elts[0];
            indA1 = cn1[e1]-1;
            indB1 = cn2[e1]-1;
            indp = indA1;
            indm = indB1;
          }
          else
          {
            indp = ind;
            indm = ind;
          }  
          ptA1[0] = xt[ind]; ptA1[1] = yt[ind]; ptA1[2] = zt[ind];
          ptB1[0] = xt[indm]; ptB1[1] = yt[indm]; ptB1[2] = zt[indm];
          ptA2[0] = xt[ind]; ptA2[1] = yt[ind]; ptA2[2] = zt[ind];
          ptB2[0] = xt[indp]; ptB2[1] = yt[indp]; ptB2[2] = zt[indp];

          alpha0 = K_COMPGEOM::getAlphaAngleBetweenBars(ptA1, ptB1, ptA2, ptB2, dirVect);
          if (alpha0 == -1000.) alpt[ind] = 0.;
          else alpt[ind] = alpha0;
        }

        RELEASESHAREDU(array,f, cn); return tpl;  
      }
    }
  }
}
