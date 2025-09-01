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

extern "C" 
{
  void k6conv2center1_(const E_Int& ni, const E_Int& nj, const E_Int& nk, 
                       const E_Int& nfld, E_Float* fieldnode, 
                       E_Float* fieldcenter);

  void k6compstructgrad_(
    const E_Int& ni, const E_Int& nj, const E_Int& nk, 
    const E_Int& nbcell, const E_Int& nbint, 
    const E_Float* xt, const E_Float* yt, const E_Float* zt, 
    const E_Float* field,
    E_Float* gradx, E_Float* grady, E_Float* gradz,
    E_Float* surf, E_Float* snorm, E_Float* centerInt,
    E_Float* vol, E_Float* fieldint);

  void k6compstructgrad1d_(
    const E_Int& ni, const E_Int& nbcell, 
    const E_Float* xt, const E_Float* yt, const E_Float* zt, 
    const E_Float* field,
    E_Float* gradx, E_Float* grady, E_Float* gradz);

  void k6compstructgrad2d_(
    const E_Int& ni, const E_Int& nj, const E_Int& nbcell, 
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    const E_Float* field,
    E_Float* surf, E_Float* nxt, E_Float* nyt, E_Float* nzt,
    E_Float* gradx, E_Float* grady, E_Float* gradz);
}

//=============================================================================
/* Calcul du gradient d un ensemble de champs definis en noeuds. 
   Le gradient est fourni aux centres des cellules */
//=============================================================================
PyObject* K_POST::computeGrad(PyObject* self,PyObject* args)
{
  PyObject* array; PyObject* varname;
  if (!PYPARSETUPLE_(args, OO_, &array, &varname)) return NULL;
  
  // Check array
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cn;
  E_Int ni, nj, nk;// number of points of array
  E_Int posx = -1; E_Int posy = -1; E_Int posz = -1;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, cn, 
                                     eltType);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "computeGrad: invalid array.");
    return NULL;
  }
  // check varname
  char* var = NULL;
  if (PyString_Check(varname)) var = PyString_AsString(varname);
#if PY_VERSION_HEX >= 0x03000000
  else if (PyUnicode_Check(varname)) var = (char*)PyUnicode_AsUTF8(varname);
#endif
  else
  {    
    PyErr_SetString(PyExc_TypeError, 
                    "computeGrad: varname must be a string.");
    RELEASESHAREDB(res,array,f,cn); return NULL;
  } 
  E_Int posv = K_ARRAY::isNamePresent(var, varString);
  if (posv == -1)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "computeGrad: variable not found in array.");
    RELEASESHAREDB(res,array,f,cn); return NULL;
  }
  posv++;
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "computeGrad: coordinates not found in array.");
    RELEASESHAREDB(res,array,f,cn); return NULL;
  }
  posx++; posy++; posz++;
  E_Int nfld = f->getNfld();
  if (nfld == 3)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "computeGrad: no field to compute.");
    RELEASESHAREDB(res,array, f, cn); return NULL;
  }
  PyObject* tpl = NULL;
  char* varStringOut;
  computeGradVarsString(var, varStringOut);

  if (res == 1) 
  {
    E_Int ni1 = K_FUNC::E_max(1,ni-1);
    E_Int nj1 = K_FUNC::E_max(1,nj-1);
    E_Int nk1 = K_FUNC::E_max(1,nk-1);    
    E_Int ncells = ni1*nj1*nk1;
    tpl = K_ARRAY::buildArray3(3, varStringOut, ni1, nj1, nk1);
    E_Float* fnp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF fp(ncells, 3, fnp, true); 
    computeGradStruct(ni, nj, nk, 
                      f->begin(posx), f->begin(posy), f->begin(posz), 
                      f->begin(posv),
                      fp.begin(1), fp.begin(2), fp.begin(3));      
  }
  else if (res == 2) 
  {
    if (strcmp(eltType, "NGON") == 0)
    {
      // Build unstructured NGON array from existing connectivity & empty fields
      E_Int nelts = cn->getNElts();
      FldArrayF* fp = new FldArrayF(nelts, 3, true); fp->setAllValuesAtNull();
      tpl = K_ARRAY::buildArray3(*fp, varStringOut, *cn, "NGON");
      delete fp; K_ARRAY::getFromArray3(tpl, fp);
      
      E_Int err = computeGradNGon(
        f->begin(posx), f->begin(posy), f->begin(posz), 
        f->begin(posv), *cn, 
        fp->begin(1), fp->begin(2), fp->begin(3));
      if (err == 1)
      {
        PyErr_SetString(PyExc_TypeError, 
                        "computeGrad: gradient can only be computed for 3D NGONs.");
        delete [] varStringOut;
        RELEASESHAREDB(res,array,f,cn); return NULL;         
      }
      RELEASESHAREDS(tpl, fp);
    }
    else  // ME
    {
      E_Int nc = cn->getNConnect();
      E_Int npts = f->getSize();
      E_Int api = f->getApi();
      E_Bool center = true;
      E_Bool copyConnect = true;
      if (nc > 1) api = 3;
      
      tpl = K_ARRAY::buildArray3(
        3, varStringOut, npts, *cn, eltType,
        center, api, copyConnect);
      FldArrayF* f2;
      K_ARRAY::getFromArray3(tpl, f2);
      E_Float* gradx = f2->begin(1);
      E_Float* grady = f2->begin(2);
      E_Float* gradz = f2->begin(3);
      computeGradUnstr(
        f->begin(posx), f->begin(posy), f->begin(posz), *cn, eltType,
        f->begin(posv), gradx, grady, gradz
      );
      RELEASESHAREDS(tpl, f2);
    }
  }
  RELEASESHAREDB(res, array, f, cn);
  delete [] varStringOut;
  return tpl;
}

//=============================================================================
/* gradx, grady, gradz must be allocated before */
//=============================================================================
E_Int K_POST::computeGradStruct(E_Int ni, E_Int nj, E_Int nk, 
                                E_Float* xt, E_Float* yt, E_Float* zt, 
                                E_Float* field,
                                E_Float* gradx, E_Float* grady, E_Float* gradz)
{
  E_Int ni1 = ni-1; E_Int nj1 = nj-1; E_Int nk1 = nk-1;
  E_Int dim = 3;
  if (ni1*nj1*nk1 == 0)
  {
    if ((ni1*nj1 == 0) && (ni1*nk1 == 0) && (nj1*nk1 == 0))
    {
      if ((ni == 1) && (nj == 1))
      {
        ni = nk;
      }
      else if ((ni == 1) && (nk == 1))
      {
        ni = nj;
      }
      dim = 1;
      nj = 2; nk = 2;
    }
    else
    {
      dim = 2;
      if (ni == 1)
      {
        ni = nj; nj = nk;
      }
      else if (nj == 1)
      {
        nj = nk;
      }
      nk = 2;
    }
  }
  ni1 = K_FUNC::E_max(1,ni-1);
  nj1 = K_FUNC::E_max(1,nj-1);
  nk1 = K_FUNC::E_max(1,nk-1);    
  E_Int ncells = ni1*nj1*nk1;
  E_Int nint = ni*nj1*nk1 + ni1*nj*nk1 + ni1*nj1*nk;
  // Construction des tableaux locaux
  if (dim == 1)
  {
    k6compstructgrad1d_(ni, ni1, xt, yt, zt, field, gradx, grady, gradz);
  }
  else if (dim == 2)
  {
    FldArrayF surf(ncells);
    FldArrayF nxt(ncells);
    FldArrayF nyt(ncells);
    FldArrayF nzt(ncells);
    k6compstructgrad2d_(ni, nj, ni1*nj1, xt, yt, zt, field, 
                        surf.begin(), nxt.begin(), nyt.begin(), nzt.begin(), 
                        gradx, grady, gradz);
  }
  else //3d
  {
    FldArrayF surf(nint, 3);
    FldArrayF snorm(nint);
    FldArrayF centerInt(nint, 3);
    FldArrayF vol(ncells);
    FldArrayF fieldint(nint);
    k6compstructgrad_(ni, nj, nk, ncells, nint, 
                      xt, yt, zt, field, gradx, grady, gradz,
                      surf.begin(), snorm.begin(), centerInt.begin(), 
                      vol.begin(), fieldint.begin());
  }
  return 1;
}

//=============================================================================
/* A partir de la chaine de variables initiale: (x,y,z,var1,var2,...)
   Cree la chaine (gradxvar1,gradyvar1,gradzvar1, gradxvar2, ....) 
   Cette routine alloue varStringOut */
//=============================================================================
void K_POST::computeGradVarsString(char* varString, char*& varStringOut)
{
  vector<char*> vars;
  K_ARRAY::extractVars(varString, vars);
  E_Int c = -1;
  E_Int varsSize = vars.size();
  E_Int sizeVarStringOut = 0;
  for (E_Int v = 0; v < varsSize; v++)
  {
    E_Int vsize = strlen(vars[v]);
    sizeVarStringOut+=vsize+6;//gradxvarString,
  }
  varStringOut = new char [3*sizeVarStringOut];

  for (E_Int v = 0; v < varsSize; v++)
  {
    char*& var0 = vars[v];
    if (strcmp(var0, "x") != 0 && 
        strcmp(var0, "y") != 0 && 
        strcmp(var0, "z") != 0)
    {
      if (c == -1)
      {
        strcpy(varStringOut, "gradx");
        c = 1;
      }
      else strcat(varStringOut, ",gradx");
      strcat(varStringOut, var0);
      strcat(varStringOut, ",grady");
      strcat(varStringOut, var0);
      strcat(varStringOut, ",gradz");
      strcat(varStringOut, var0);
    }
  } 
  for (E_Int v = 0; v < varsSize; v++) delete [] vars[v];
}

//==============================================================================
E_Int K_POST::computeGradNGon(E_Float* xt, E_Float* yt, E_Float* zt, 
                              E_Float* fp, FldArrayI& cn,
                              E_Float* gradx, E_Float* grady, E_Float* gradz)
{
  // Acces non universel sur le ptrs
  E_Int* ngon = cn.getNGon();
  E_Int* nface = cn.getNFace();
  E_Int* indPG = cn.getIndPG();
  E_Int* indPH = cn.getIndPH();
  // Acces universel nbres d'elements et de faces
  E_Int nelts = cn.getNElts();
  E_Int nfaces = cn.getNFaces();
  E_Int ierr = 0; // error index

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
  vector< vector<E_Int> > cnEV(nelts);
  K_CONNECT::connectNG2EV(cn, cnEV); //deja calculee dans NGONVol
  // Tableau de la dimension des elements
  FldArrayI dimElt(nelts); 
  K_CONNECT::getDimElts(cn, dimElt);

#pragma omp parallel
  {
    E_Float fpmeanface, invvol;
    E_Int dim, ind, noface, indnode, nbFaces, nbNodes, nbNodesPerFace;
    E_Float xbe, ybe, zbe; // coordonnees du barycentre d un element
    E_Float xbf, ybf, zbf; // coordonnees du barycentre d une face
    E_Float sens, sx, sy, sz;

    // parcours des elements
#pragma omp for
    for (E_Int et = 0; et < nelts; et++)
    { 
      // An error occurred - skip the rest of the calculations
      if (ierr == 1) continue;
      
      dim = dimElt[et]; // dimension de l'element
      if (dim == 3) 
      {
        invvol = 1./volp[et];
        gradx[et] = 0.; grady[et] = 0.; gradz[et] = 0.;
        
        // calcul du barycentre be (xbe, ybe, zbe) de l'element
        const vector<E_Int>& vertices = cnEV[et]; // sommets associes a l'element et
        nbNodes = vertices.size();
        xbe = 0.; ybe = 0.; zbe = 0.;
        for (E_Int n = 0; n < nbNodes; n++)
        {
          ind = vertices[n]-1;
          xbe += xt[ind]; ybe += yt[ind]; zbe += zt[ind];
        }
        xbe = xbe/nbNodes; ybe = ybe/nbNodes; zbe = zbe/nbNodes;

        // Acces universel element et
        E_Int* elt = cn.getElt(et, nbFaces, nface, indPH);
        // parcours des faces de l element et
        for (E_Int fa = 0; fa < nbFaces; fa++)
        {
          // Acces universel face noface
          noface = elt[fa]-1;
          E_Int* face = cn.getFace(noface, nbNodesPerFace, ngon, indPG);
          //valeur moyenne de fp pour la face
          fpmeanface = 0.;
          // calcul du barycentre bf (xbf, ybf, zbf) de la face
          xbf = 0.; ybf = 0.; zbf = 0.;
          for (E_Int n = 0; n < nbNodesPerFace; n++)
          {
            indnode = face[n]-1; // indice du point n
            xbf += xt[indnode]; ybf += yt[indnode]; zbf += zt[indnode];
            fpmeanface += fp[indnode];
          }
          xbf = xbf/nbNodesPerFace; ybf = ybf/nbNodesPerFace; zbf = zbf/nbNodesPerFace;            
          fpmeanface = fpmeanface/nbNodesPerFace;           
          // bilan
          // verification du sens de la normale. Celle-ci doit etre exterieure
          sx = sxp[noface]; sy = syp[noface]; sz = szp[noface];
          sens = (xbe-xbf)*sx + (ybe-ybf)*sy + (zbe-zbf)*sz;
          if (sens > 0.) {sx=-sx; sy=-sy; sz=-sz;}
          gradx[et] += fpmeanface*sx;
          grady[et] += fpmeanface*sy;
          gradz[et] += fpmeanface*sz;
        }
        gradx[et] *= invvol;
        grady[et] *= invvol;
        gradz[et] *= invvol;
      }
      else
      {
        printf("computeGrad: not valid for " SF_D_ "D NGONs\n", dim);
        ierr = 1;
      }
    }
  }
  delete [] volp;
  delete [] sxp; 
  delete [] syp;
  delete [] szp;
  delete [] snp;
  return ierr;
}
