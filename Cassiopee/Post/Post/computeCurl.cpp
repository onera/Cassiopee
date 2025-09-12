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
# include <string.h>
# include "post.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
/* Calcul du rotationnel d'un champ defini par un vecteur (u,v,w) en noeuds. 
   Le rotationnel est fourni aux centres des cellules */
//=============================================================================
PyObject* K_POST::computeCurl(PyObject* self, PyObject* args)
{
  PyObject* array; PyObject* vars0;
  if (!PYPARSETUPLE_(args, OO_, &array, &vars0)) return NULL;
  
  //extraction des variables constituant le vecteur dont le rot est calcule
  vector<char*> vars;
  if (PyList_Check(vars0) == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeCurl: a list of 3 variables for curl computation must be defined.");
    return NULL; 
  }
  if (PyList_Size(vars0) != 3)
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeCurl: 3 variables must be defined to extract the curl.");
    return NULL;
  }
  for (Py_ssize_t i = 0; i < PyList_Size(vars0); i++)
  {
    PyObject* tpl0 = PyList_GetItem(vars0, i);
    if (PyString_Check(tpl0))
    {
      char* str = PyString_AsString(tpl0);
      vars.push_back(str);
    }
#if PY_VERSION_HEX >= 0x03000000
    else if (PyUnicode_Check(tpl0))
    {
      char* str = (char*)PyUnicode_AsUTF8(tpl0);
      vars.push_back(str);  
    }
#endif
    else  
    {
      PyErr_SetString(PyExc_TypeError,
                      "computeCurl: varname must be a string.");
      return NULL;
    }
  }
  // Check array
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cn;
  E_Int ni, nj, nk; // number of points of array
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, cn, eltType);
  
  if (res != 1 && res != 2) 
  {
    PyErr_SetString(PyExc_TypeError, 
                    "computeCurl: invalid array.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }

  E_Int posu = K_ARRAY::isNamePresent(vars[0], varString);
  E_Int posv = K_ARRAY::isNamePresent(vars[1], varString);
  E_Int posw = K_ARRAY::isNamePresent(vars[2], varString);
  if (posu == -1 || posv == -1 || posw == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeCurl: one variable was not found in array.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }
  posu++; posv++; posw++;
  
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "computeCurl: coordinates not found in array.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }
  
  posx++; posy++; posz++;
  E_Int nfld = f->getNfld();
  if (nfld < 6)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "computeCurl: no field to compute.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }

  // calcul du rotationnel 
  PyObject* tpl;
  char* varStringOut = new char[15];
  strcpy(varStringOut, "rotx,roty,rotz"); 
  if (res == 1)
  {
    E_Int ni1 = K_FUNC::E_max(1, ni-1);
    E_Int nj1 = K_FUNC::E_max(1, nj-1);
    E_Int nk1 = K_FUNC::E_max(1, nk-1);
    E_Int ncells = ni1*nj1*nk1;
    tpl = K_ARRAY::buildArray(3, varStringOut, ni1, nj1, nk1);
    E_Float* fnp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF fp(ncells, 3, fnp, true); 
    computeCurlStruct(ni, nj, nk, 
                      f->begin(posx), f->begin(posy), f->begin(posz),
                      f->begin(posu), f->begin(posv), f->begin(posw),
                      fp.begin(1), fp.begin(2), fp.begin(3));
  }
  else // non structure
  {
    // if (strcmp(eltType, "NGON") == 0)
    // {
    //   E_Int api = f->getApi();
    //   E_Int npts = f->getSize();
    //   E_Int sizeFN = cn->getSizeNGon(); //  taille de la connectivite Face/Noeuds
    //   E_Int nelts = cn->getNElts(); // nombre total d elements
    //   E_Int nfaces = cn->getNFaces();
    //   tpl = K_ARRAY::buildArray3(3, varStringOut, npts, nelts, nfaces,
    //                              eltType, true, api);
    //   FldArrayF* fnp; FldArrayI* cnnp;
    //   K_ARRAY::getFromArray3(tpl, fnp, cnnp);
    //   FldArrayI cnn(cn->getSize(), 1, cnnp, true); cnn = *cn;
    //   E_Int ierr = computeCurlNGon(
    //     f->begin(posx), f->begin(posy), f->begin(posz),
    //     f->begin(posu), f->begin(posv), f->begin(posw), *cn, 
    //     fp.begin(1), fp.begin(2), fp.begin(3)
    //   );
    //   RELEASESHAREDU(array, f, cn);

    //   if (ierr == 1)
    //   {
    //     PyErr_SetString(PyExc_TypeError, 
    //                     "computeCurl: curl can only be computed for 3D NGONs.");
    //     RELEASESHAREDB(res, array, f, cn); return NULL;         
    //   }      
    // }
    
    if (strcmp(eltType, "NGON") == 0)
    {
      E_Int npts = f->getSize();
      E_Int* cnp = cn->begin();
      //E_Int nfaces = cnp[0]; // nombre total de faces
      E_Int sizeFN = cnp[1]; //  taille de la connectivite Face/Noeuds
      E_Int nelts = cnp[sizeFN+2];  // nombre total d elements
      E_Int csize = cn->getSize();
      tpl = K_ARRAY::buildArray(3, varStringOut, npts, nelts, -1, eltType, true, csize);
      E_Float* fnp = K_ARRAY::getFieldPtr(tpl);
      FldArrayF fp(nelts, 3, fnp, true);
      E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
      FldArrayI cnn(cn->getSize(), 1, cnnp, true); cnn = *cn;
      E_Int err = computeCurlNGon(
        f->begin(posx), f->begin(posy), f->begin(posz),
        f->begin(posu), f->begin(posv), f->begin(posw), *cn, 
        fp.begin(1), fp.begin(2), fp.begin(3));    

      if (err == 1)
      {
        PyErr_SetString(PyExc_TypeError, 
                        "computeCurl: curl can only be computed for 3D NGONs.");
        RELEASESHAREDB(res, array, f, cn); return NULL;         
      }      
    }
    else if (strcmp(eltType, "TRI") == 0 ||
             strcmp(eltType, "QUAD")  == 0 ||
             strcmp(eltType, "TETRA") == 0 ||
             strcmp(eltType, "HEXA") == 0 ||
             strcmp(eltType, "PENTA") == 0) 
    {
      E_Int npts = f->getSize();
      tpl = K_ARRAY::buildArray(3, varStringOut, npts, cn->getSize(),-1,eltType,true,cn->getSize()*cn->getNfld());
      E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
      K_KCORE::memcpy__(cnnp, cn->begin(), cn->getSize()*cn->getNfld());
      E_Float* fnp = K_ARRAY::getFieldPtr(tpl);
      E_Int nelts = cn->getSize();    
      FldArrayF fp(nelts, 3, fnp, true);
      
      // calcul du rotationnel aux centres des elements
      computeCurlUnstruct(
        *cn, eltType,
        f->begin(posx), f->begin(posy), f->begin(posz),
        f->begin(posu), f->begin(posv), f->begin(posw),
        fp.begin(1), fp.begin(2), fp.begin(3)
      );              
    }
    else
    {
      PyErr_SetString(PyExc_TypeError,
                      "computeCurl: not a valid element type.");
      RELEASESHAREDU(array, f, cn); return NULL;
    }    
  }
  
  RELEASESHAREDB(res, array, f, cn);
  delete [] varStringOut;
  return tpl;
}
//==============================================================================
E_Int K_POST::computeCurlStruct(
  const E_Int ni, const E_Int nj, const E_Int nk, 
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  const E_Float* ux, const E_Float* uy, const E_Float* uz,
  E_Float* rotx, E_Float* roty, E_Float* rotz
)
{
  if (ni*nj == 1 || ni*nk == 1 || nj*nk == 1) return -1;
  if (ni == 1 || nj == 1 || nk == 1)
  {
    compCurlStruct2D(ni, nj, nk, xt, yt, zt, ux, uy, uz, rotx, roty, rotz);
  }
  else
  {
    compCurlStruct3D(ni, nj, nk, xt, yt, zt, ux, uy, uz, rotx, roty, rotz);
  }
  return 1;
}
//=============================================================================
E_Int K_POST::computeCurlUnstruct(
  FldArrayI& cn, const char* eltType,
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  const E_Float* ux, const E_Float* uy, const E_Float* uz,
  E_Float* rotx, E_Float* roty, E_Float* rotz
)
{
  // Get ME mesh dimensionality from the first element type
  E_Int dim = 3;
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);
  if (strcmp(eltTypes[0], "BAR") == 0) dim = 1;
  else if (strcmp(eltTypes[0], "TRI") == 0 or
           strcmp(eltTypes[0], "QUAD") == 0) dim = 2;
  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];

  if (dim == 2)
  {
    compCurlUnstruct2D(cn, eltType, xt, yt, zt, ux, uy, uz, rotx, roty, rotz);
  }
  else if (dim == 3)
  {
    compCurlUnstruct3D(cn, eltType, xt, yt, zt, ux, uy, uz, rotx, roty, rotz);
  }
  else return -1;
  return 1;
}
//==============================================================================
E_Int K_POST::computeCurlNGon(
  const E_Float* xt, const E_Float* yt, const E_Float* zt, 
  const E_Float* fxp, const E_Float* fyp, const E_Float* fzp, FldArrayI& cn,
  E_Float* curlx, E_Float* curly, E_Float* curlz
)
{
  // Donnees liees a la connectivite
  E_Int nfaces = cn.getNFaces(); // nombre total de faces
  E_Int nelts = cn.getNElts();  // nombre total d elements
  E_Int* ngon = cn.getNGon(); E_Int* indPG = cn.getIndPG();
  E_Int* nface = cn.getNFace(); E_Int* indPH = cn.getIndPH();

  // calcul de la metrique
  E_Float* sxp = new E_Float [3*nfaces];
  E_Float* syp = new E_Float [3*nfaces];
  E_Float* szp = new E_Float [3*nfaces];
  E_Float* snp = new E_Float [nfaces];
  FldArrayI* cFE = new FldArrayI();
  K_CONNECT::connectNG2FE(cn, *cFE);
  K_METRIC::compNGonFacesSurf(xt, yt, zt, cn, sxp, syp, szp, snp, cFE);
  delete cFE;
  E_Float* volp = new E_Float [nelts];
  K_METRIC::compNGonVol(xt, yt, zt, cn, volp); 
  // Connectivite Element/Noeuds
  vector<vector<E_Int> > cnEV(nelts);
  K_CONNECT::connectNG2EV(cn, cnEV); //deja calculee dans NGONVol

  FldArrayI dimElt(nelts); // tableau de la dimension des elements
  K_CONNECT::getDimElts(cn, dimElt);
  if (dimElt[0] < 3)
  {
    printf("computeCurl: not valid for " SF_D_ "D NGONs\n", dimElt[0]);
    delete [] volp;
    delete [] sxp; 
    delete [] syp;
    delete [] szp;
    delete [] snp;
    return 1;
  }

  #pragma omp parallel
  {
    E_Int ind, noface, indnode, nbFaces, nbNodes, nbNodesPerFace;
    E_Float fxpmeanface, fypmeanface, fzpmeanface, invvol;
    E_Float xbe, ybe, zbe; // coordonnees du barycentre d un element
    E_Float xbf, ybf, zbf; // coordonnees du barycentre d une face
    E_Float sens, sx, sy, sz;
    vector<E_Int> vertices; // sommets associes a l'element

    // parcours des elements
    #pragma omp for
    for (E_Int et = 0; et < nelts; et++)
    { 
      invvol = -1./volp[et];
      curlx[et] = 0.; curly[et] = 0.; curlz[et] = 0.;
      
      // calcul du barycentre be (xbe, ybe, zbe) de l'element
      vertices = cnEV[et];
      nbNodes = vertices.size();
      xbe = 0.; ybe = 0.; zbe = 0.;
      for (E_Int n = 0; n < nbNodes; n++)
      {
        ind = vertices[n]-1;
        xbe += xt[ind]; ybe += yt[ind]; zbe += zt[ind];
      }
      xbe = xbe/nbNodes; ybe = ybe/nbNodes; zbe = zbe/nbNodes;

      // parcours des faces de l element et
      E_Int* elt = cn.getElt(et, nbFaces, nface, indPH);
      for (E_Int fa = 0; fa < nbFaces; fa++)
      {
        noface = elt[fa]-1;
        E_Int* face = cn.getFace(noface, nbNodesPerFace, ngon, indPG);
        // valeur moyenne de fp pour la face
        fxpmeanface = 0.; fypmeanface = 0.; fzpmeanface = 0.;
        // calcul du barycentre bf (xbf, ybf, zbf) de la face
        xbf = 0.; ybf = 0.; zbf = 0.;
        for (E_Int n = 0; n < nbNodesPerFace; n++)
        {
          indnode = face[n]-1;
          xbf += xt[indnode]; ybf += yt[indnode]; zbf += zt[indnode];
          fxpmeanface += fxp[indnode];
          fypmeanface += fyp[indnode];
          fzpmeanface += fzp[indnode];
        }
        xbf = xbf/nbNodesPerFace; ybf = ybf/nbNodesPerFace; zbf = zbf/nbNodesPerFace;            
        fxpmeanface = fxpmeanface/nbNodesPerFace;           
        fypmeanface = fypmeanface/nbNodesPerFace;           
        fzpmeanface = fzpmeanface/nbNodesPerFace;           

        // bilan
        // verification du sens de la normale. Celle-ci doit etre exterieure
        sx = sxp[noface]; sy = syp[noface]; sz = szp[noface];
        sens = (xbe-xbf)*sx + (ybe-ybf)*sy + (zbe-zbf)*sz;
        if (sens > 0.) {sx=-sx; sy=-sy; sz=-sz;}
        curlx[et] += fypmeanface*sz - fzpmeanface*sy;
        curly[et] += fzpmeanface*sx - fxpmeanface*sz;
        curlz[et] += fxpmeanface*sy - fypmeanface*sx;
      }
      curlx[et] *= invvol;
      curly[et] *= invvol;
      curlz[et] *= invvol;
    }
  }

  delete [] volp;
  delete [] sxp; 
  delete [] syp;
  delete [] szp;
  delete [] snp;
  return 0;
}
