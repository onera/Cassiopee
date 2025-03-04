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

  void k6compstructdiv_(
    const E_Int& ni, const E_Int& nj, const E_Int& nk,
    const E_Int& nbcell, const E_Int& nbint,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    const E_Float* fieldX, const E_Float* fieldY, const E_Float* fieldZ,
    E_Float* div,
    E_Float* surf, E_Float* snorm, E_Float* centerInt,
    E_Float* vol, E_Float* fldintX, E_Float* fldintY, E_Float* fldintZ);

  void k6compstructdiv2d_(
    const E_Int& ni, const E_Int& nj, const E_Int& nbcell,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    const E_Float* fieldX, const E_Float* fieldY, const E_Float* fieldZ,
    E_Float* surf, E_Float* nxt, E_Float* nyt, E_Float* nzt,
    E_Float* div);

  void k6compunstrdiv_(E_Int& dim, E_Int& npts, E_Int& nelts,
                       E_Int& nedges, E_Int& nnodes,
                       E_Int* cn, E_Float* xt, E_Float* yt, E_Float* zt,
                       E_Float* fieldX, E_Float* fieldY, E_Float* fieldZ,
                       E_Float* fieldfX, E_Float* fieldfY, E_Float* fieldfZ,
                       E_Float* snx, E_Float* sny, E_Float* snz,
                       E_Float* surf, E_Float* vol,
                       E_Float* xint, E_Float* yint, E_Float* zint,
                       E_Float* div);

  void k6compunstrdiv2d_(E_Int& npts, E_Int& nelts, E_Int& nnodes,
                         E_Int* cn, E_Float* xt, E_Float* yt, E_Float* zt,
                         E_Float* fieldX, E_Float* fieldY, E_Float* fieldZ,
                         E_Float* snx, E_Float* sny, E_Float* snz,
                         E_Float* surf, E_Float* div);
}

//=============================================================================
/* Compute the divergence of a set of vector fields given at mesh nodes.
   The divergence is given on cell centers. */
//=============================================================================
PyObject* K_POST::computeDiv(PyObject* self,PyObject* args)
{
  PyObject* array; PyObject* vars0;
  if (!PyArg_ParseTuple(args, "OO", &array, &vars0)) return NULL;

  //extract variables constituting components of the vector whose div is calculated
  vector<char*> vars;
  if (PyList_Check(vars0) == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeDiv: a list of 3 variables for div computation must be defined.");
    return NULL;
  }
  if (PyList_Size(vars0) != 3)
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeDiv: 3 variables must be defined to extract the div.");
    return NULL;
  }
  for (int i = 0; i < PyList_Size(vars0); i++)
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
                      "computeDiv: varname must be a string.");
      return NULL;
    }
  }

  // Check array
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cn;
  E_Int ni, nj, nk;// number of points of array
  E_Int posx = -1; E_Int posy = -1; E_Int posz = -1;
  E_Int res = K_ARRAY::getFromArray(array, varString, f, ni, nj, nk, cn,
                                    eltType, true);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeDiv: invalid array.");
    return NULL;
  }
  E_Int posu = K_ARRAY::isNamePresent(vars[0], varString);
  E_Int posv = K_ARRAY::isNamePresent(vars[1], varString);
  E_Int posw = K_ARRAY::isNamePresent(vars[2], varString);
  if (posu == -1 || posv == -1 || posw == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeDiv: at least one variable was not found in array.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }
  posu++; posv++; posw++;

  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeDiv: coordinates not found in array.");
    RELEASESHAREDB(res,array,f,cn); return NULL;
  }
  posx++; posy++; posz++;
  E_Int nfld = f->getNfld();
  if (nfld < 6)
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeDiv: no field to compute.");
    RELEASESHAREDB(res,array, f, cn); return NULL;
  }

  if (strncmp(vars[0], vars[1], strlen(vars[0])-1) != 0 ||
      strncmp(vars[1], vars[2], strlen(vars[1])-1) != 0 )
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeDiv: invalid names for vector component fields.");
    RELEASESHAREDB(res,array,f,cn); return NULL;
  }
  char* sv0 = vars[0]; char* sv1 = vars[1]; char* sv2 = vars[2];
  char s0 = sv0[strlen(sv0)-1];
  char s1 = sv1[strlen(sv1)-1];
  char s2 = sv2[strlen(sv2)-1];
  if (s0 != 'X' || s1 != 'Y' || s2 != 'Z')
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeDiv: error with the order of given scalar fields.");
    RELEASESHAREDB(res,array,f,cn); return NULL;
  }
  char* varStringOut = new char [strlen(vars[0])+3]; // +3 for "div", -1 for trailing 'X'
  strcpy(varStringOut, "div");
  char* pt = varStringOut;
  char* v = vars[0];
  pt += 3;
  for (size_t j = 0; j < strlen(v)-1; j++)
  {
    *pt = v[j]; pt++;
  }
  *pt = '\0';

  PyObject* tpl = NULL;
  if (res == 1)
  {
    E_Int ni1 = K_FUNC::E_max(1,ni-1);
    E_Int nj1 = K_FUNC::E_max(1,nj-1);
    E_Int nk1 = K_FUNC::E_max(1,nk-1);
    E_Int ncells = ni1*nj1*nk1;
    tpl = K_ARRAY::buildArray(1, varStringOut, ni1, nj1, nk1);
    E_Float* fnp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF fp(ncells, 1, fnp, true); fp.setAllValuesAtNull();
    E_Int err = computeDivStruct(ni, nj, nk,
                                 f->begin(posx), f->begin(posy), f->begin(posz),
                                 f->begin(posu), f->begin(posv), f->begin(posw),
                                 fp.begin());
    if (err == -1)
    {
      PyErr_SetString(PyExc_TypeError,
                      "computeDiv: Not valid for 1D.");
      RELEASESHAREDB(res,array,f,cn); return NULL;
    }
  }
  else if (res == 2)
  {
    if (strcmp(eltType, "NGON") == 0)
    {
      E_Int npts = f->getSize();
      E_Int* cnp = cn->begin();
      //E_Int nfaces = cnp[0]; // nombre total de faces
      E_Int sizeFN = cnp[1]; //  taille de la connectivite Face/Noeuds
      E_Int nelts = cnp[sizeFN+2];  // nombre total d elements
      E_Int csize = cn->getSize();
      tpl = K_ARRAY::buildArray(1, varStringOut, npts, nelts, -1, eltType, true, csize);
      E_Float* fnp = K_ARRAY::getFieldPtr(tpl);
      FldArrayF fp(nelts, 1, fnp, true);
      E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
      FldArrayI cnn(cn->getSize(), 1, cnnp, true); cnn = *cn;
      E_Int err = computeDivNGon(f->begin(posx), f->begin(posy), f->begin(posz),
                                 f->begin(posu), f->begin(posv), f->begin(posw), *cn,
                                 fp.begin());
      if (err == 1)
      {
        PyErr_SetString(PyExc_TypeError,
                        "computeDiv: divergence can only be computed for 3D NGONs.");
        RELEASESHAREDB(res,array,f,cn); return NULL;
      }
    }
    else
    {
      if (strcmp(eltType, "TRI") != 0 &&
          strcmp(eltType, "QUAD") != 0 &&
          strcmp(eltType, "TETRA") != 0 &&
          strcmp(eltType, "HEXA") != 0 &&
          strcmp(eltType, "PENTA") != 0 )
      {
        PyErr_SetString(PyExc_TypeError,
                        "computeDiv: not a valid element type.");
        RELEASESHAREDU(array,f,cn); return NULL;
      }
      E_Int npts = f->getSize();
      tpl = K_ARRAY::buildArray(1, varStringOut, npts, cn->getSize(), -1, eltType, true,
                                cn->getSize()*cn->getNfld());
      E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
      K_KCORE::memcpy__(cnnp, cn->begin(), cn->getSize()*cn->getNfld());
      E_Float* fnp = K_ARRAY::getFieldPtr(tpl);
      E_Int nelts = cn->getSize();
      FldArrayF fp(nelts, 1, fnp, true);
      computeDivNS(eltType, npts, *cn,
                   f->begin(posx), f->begin(posy), f->begin(posz),
                   f->begin(posu), f->begin(posv), f->begin(posw),
                   fp.begin());
    }
  }
  RELEASESHAREDB(res,array,f,cn);
  delete [] varStringOut;
  return tpl;
}

//=============================================================================
/* div must be allocated before */
//=============================================================================
E_Int K_POST::computeDivStruct(E_Int ni, E_Int nj, E_Int nk,
                               E_Float* xt, E_Float* yt, E_Float* zt,
                               E_Float* fieldX, E_Float* fieldY, E_Float* fieldZ,
                               E_Float* div)
{
  if ((ni == 1 && nj == 1) || (ni == 1 && nk == 1) || (nj == 1 && nk == 1))
    return -1;

  E_Int dim = 3;
  if (ni == 1)
  {
    ni = nj; nj = nk; nk = 2; dim = 2;
  }
  else if (nj == 1)
  {
    nj = nk; nk = 2; dim = 2;
  }
  else if (nk == 1)
  {
    nk = 2; dim = 2;
  }
  E_Int ni1 = K_FUNC::E_max(1,ni-1);
  E_Int nj1 = K_FUNC::E_max(1,nj-1);
  E_Int nk1 = K_FUNC::E_max(1,nk-1);
  E_Int ncells = ni1*nj1*nk1;
  E_Int nint = ni*nj1*nk1 + ni1*nj*nk1 + ni1*nj1*nk;
  // Construction des tableaux locaux
  if (dim == 2)
  {
    FldArrayF surf(ncells);
    FldArrayF nxt(ncells);
    FldArrayF nyt(ncells);
    FldArrayF nzt(ncells);
    k6compstructdiv2d_(ni, nj, ni1*nj1, xt, yt, zt, fieldX, fieldY, fieldZ,
                       surf.begin(), nxt.begin(), nyt.begin(), nzt.begin(),
                       div);
  }
  else //3d
  {
    FldArrayF surf(nint, 3);
    FldArrayF snorm(nint);
    FldArrayF centerInt(nint, 3);
    FldArrayF vol(ncells);
    FldArrayF fieldintX(nint);
    FldArrayF fieldintY(nint);
    FldArrayF fieldintZ(nint);
    k6compstructdiv_(ni, nj, nk, ncells, nint,
                     xt, yt, zt, fieldX, fieldY, fieldZ, div,
                     surf.begin(), snorm.begin(), centerInt.begin(),
                     vol.begin(),
                     fieldintX.begin(), fieldintY.begin(), fieldintZ.begin());
  }
  return 1;
}
//=============================================================================
E_Int K_POST::computeDivNS(char* eltType, E_Int npts, FldArrayI& cn,
                           E_Float* xt, E_Float* yt, E_Float* zt,
                           E_Float* fieldX, E_Float* fieldY, E_Float* fieldZ,
                           E_Float* div)
{
  E_Int nelts = cn.getSize();
  E_Int nnodes; //nb de noeuds par elts
  E_Int nedges; //nb de facettes par elts
  E_Int dim = 3;
  if (strcmp(eltType, "TRI") == 0)
  {
    nnodes = 3; nedges = 3; dim = 2;
  }
  else if (strcmp(eltType, "QUAD") == 0)
  {
    nnodes = 4; nedges = 4; dim = 2;
  }
  else if (strcmp(eltType, "TETRA") == 0)
  {
    nedges = 4; nnodes = 4;
  }
  else if (strcmp( eltType, "HEXA") == 0)
  {
    nedges = 6; nnodes = 8;
  }
  else if (strcmp(eltType, "PENTA") == 0)
  {
    nedges = 5; nnodes = 6;
  }
  else return -1;

  if (dim == 3)
  {
    //tmp tabs
    FldArrayF fieldfx(nelts,nedges);
    FldArrayF fieldfy(nelts,nedges);
    FldArrayF fieldfz(nelts,nedges);
    FldArrayF snx(nelts,nedges);
    FldArrayF sny(nelts,nedges);
    FldArrayF snz(nelts,nedges);
    FldArrayF surf(nelts,nedges);
    FldArrayF vol(nelts);
    FldArrayF xint(nelts, nedges);
    FldArrayF yint(nelts, nedges);
    FldArrayF zint(nelts, nedges);
    k6compunstrdiv_(dim, npts, nelts, nedges, nnodes, cn.begin(),
                    xt, yt, zt, fieldX, fieldY, fieldZ,
                    fieldfx.begin(), fieldfy.begin(), fieldfz.begin(),
                    snx.begin(), sny.begin(), snz.begin(),
                    surf.begin(), vol.begin(),
                    xint.begin(), yint.begin(), zint.begin(),
                    div);
  }
  else // dim == 2
  {
    FldArrayF snx(nelts,1);
    FldArrayF sny(nelts,1);
    FldArrayF snz(nelts,1);
    FldArrayF surf(nelts,1);
    k6compunstrdiv2d_(npts, nelts, nnodes, cn.begin(),
                      xt, yt, zt, fieldX, fieldY, fieldZ,
                      snx.begin(1), sny.begin(1), snz.begin(1),
                      surf.begin(1),
                      div);
  }
  return 1;
}

//==============================================================================
E_Int K_POST::computeDivNGon(E_Float* xt, E_Float* yt, E_Float* zt,
                             E_Float* fpx, E_Float* fpy, E_Float* fpz, FldArrayI& cn,
                             E_Float* div)
{
  E_Int* cnp = cn.begin();
  // Donnees liees a la connectivite
  E_Int nfaces = cnp[0]; // nombre total de faces
  E_Int sizeFN = cnp[1]; //  taille de la connectivite Face/Noeuds
  E_Int nelts = cnp[sizeFN+2];  // nombre total d elements
  E_Int* cEFp = cnp+4+sizeFN;// debut connectivite Elmt/Faces

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
  K_METRIC::CompNGonVol(xt, yt, zt, cn, volp);
  // Connectivite Element/Noeuds
  vector< vector<E_Int> > cnEV(nelts);
  K_CONNECT::connectNG2EV(cn, cnEV); //deja calculee dans NGONVol

  // sommets associes a l'element
  vector<E_Int> vertices;

  FldArrayI posFace(nfaces); // tableau de position des faces dans la connectivite
  K_CONNECT::getPosFaces(cn, posFace);

  E_Float fpxmeanface, fpymeanface, fpzmeanface, invvol;
  E_Int dim, ind, noface, indnode, nbFaces, nbNodes, nbNodesPerFace, pos;
  E_Float xbe, ybe, zbe; // coordonnees du barycentre d un element
  E_Float xbf, ybf, zbf; // coordonnees du barycentre d une face
  E_Float sens, sx, sy, sz;
  FldArrayI dimElt(nelts); // tableau de la dimension des elements
  K_CONNECT::getDimElts(cn, dimElt);
  // parcours des elements
  for (E_Int et = 0; et < nelts; et++)
  {
    dim = dimElt[et]; // dimension de l'element
    switch (dim)
    {
      case 1: // NGon 1D
        printf("computeDiv: not valid for 1D NGONs\n");
        delete [] volp;
        delete [] sxp;
        delete [] syp;
        delete [] szp;
        delete [] snp;
        return 1;

      case 2: // NGon 2D
        printf("computeDiv: not valid for 2D NGONs\n");
        delete [] volp;
        delete [] sxp;
        delete [] syp;
        delete [] szp;
        delete [] snp;
        return 1;

      case 3:
        invvol = 1./volp[et];
        div[et] = 0.;

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
        nbFaces = cEFp[0];
        for (E_Int fa = 0; fa < nbFaces; fa++)
        {
          noface = cEFp[fa+1]-1;
          pos = posFace[noface];
          nbNodesPerFace = cnp[pos]; pos++;
          //valeur moyenne de fpx,fpy,fpz pour la face
          fpxmeanface = 0.; fpymeanface = 0.; fpzmeanface = 0.;
          // calcul du barycentre bf (xbf, ybf, zbf) de la face
          xbf = 0.; ybf = 0.; zbf = 0.;
          for (E_Int n = 0; n < nbNodesPerFace; n++)
          {
            indnode = cnp[pos+n]-1;
            xbf += xt[indnode]; ybf += yt[indnode]; zbf += zt[indnode];
            fpxmeanface += fpx[indnode];
            fpymeanface += fpy[indnode];
            fpzmeanface += fpz[indnode];
          }
          xbf = xbf/nbNodesPerFace; ybf = ybf/nbNodesPerFace; zbf = zbf/nbNodesPerFace;
          fpxmeanface /= nbNodesPerFace;
          fpymeanface /= nbNodesPerFace;
          fpzmeanface /= nbNodesPerFace;
          // bilan
          // verification du sens de la normale. Celle-ci doit etre exterieure
          sx = sxp[noface]; sy = syp[noface]; sz = szp[noface];
          sens = (xbe-xbf)*sx + (ybe-ybf)*sy + (zbe-zbf)*sz;
          if (sens > 0.) {sx=-sx; sy=-sy; sz=-sz;}
          div[et] += fpxmeanface*sx + fpymeanface*sy + fpzmeanface*sz;
        }
        cEFp += nbFaces+1;
        div[et] *= invvol;
        break;

      default:
        delete [] volp;
        delete [] sxp;
        delete [] syp;
        delete [] szp;
        delete [] snp;
        return 1;
    }
  }
  delete [] volp;
  delete [] sxp;
  delete [] syp;
  delete [] szp;
  delete [] snp;
  return 0;
}
