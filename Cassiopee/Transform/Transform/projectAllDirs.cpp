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

# include "transform.h"
# include "Nuga/include/BbTree.h"

using namespace K_FLD;
using namespace std;

// ============================================================================
/* Projette une liste de zones 1D ou 2D sur une surface definie par array2 
   (TRI) suivant une direction donnee pour chaque point de array1 */
// ============================================================================
PyObject* K_TRANSFORM::projectAllDirs(PyObject* self, PyObject* args)
{
  PyObject* arrays; PyObject* surfArrays;
  PyObject* varsO;
  E_Int oriented;
  if (!PYPARSETUPLE_(args, OOO_ I_,
                    &arrays, &surfArrays, &varsO, &oriented))
  {
      return NULL;
  }
  // Check normal components
  if (PyList_Check(varsO) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "projectAllDirs: vars must be a list.");
    return NULL;
  }
  E_Int nvars = PyList_Size(varsO);
  if (nvars != 3)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "projectAllDirs: vars must be a 3-component vector.");
    return NULL;
  }
 
  // Extract infos from zones to be projected
  vector<E_Int> resl; vector<char*> varStringP;
  vector<FldArrayF*> fieldsP;
  vector<void*> a2; //ni,nj,nk ou cnt en NS
  vector<void*> a3; //eltType en NS
  vector<void*> a4;
  vector<PyObject*> objsP;
  E_Bool skipNoCoord = true;  E_Bool skipStructured = false;
  E_Bool skipUnstructured = false;  E_Bool skipDiffVars = true;
  E_Int isOk = K_ARRAY::getFromArrays(
    arrays, resl, varStringP, fieldsP, a2, a3, a4, objsP,  
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  E_Int nprojectedZones = fieldsP.size();
  if (isOk == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "projectAllDirs: invalid list of arrays.");
    for (E_Int no = 0; no < nprojectedZones; no++)
      RELEASESHAREDA(resl[no],objsP[no],fieldsP[no],a2[no],a3[no],a4[no]);   
    return NULL;
  }
  if (nprojectedZones == 0 )
  {
    PyErr_SetString(PyExc_TypeError,
                    "projectAllDirs: no valid projected zone found.");
    for (E_Int no = 0; no < nprojectedZones; no++)
      RELEASESHAREDA(resl[no],objsP[no],fieldsP[no],a2[no],a3[no],a4[no]);   
    return NULL;
  }

  E_Int posx1, posy1, posz1;
  vector<E_Int> posxp; vector<E_Int> posyp; vector<E_Int> poszp;
  for (E_Int no = 0; no < nprojectedZones; no++)
  {
    posx1 = K_ARRAY::isCoordinateXPresent(varStringP[no]); posx1++;
    posy1 = K_ARRAY::isCoordinateYPresent(varStringP[no]); posy1++;
    posz1 = K_ARRAY::isCoordinateZPresent(varStringP[no]); posz1++;
    posxp.push_back(posx1); posyp.push_back(posy1); poszp.push_back(posz1); 
  }

  // Normals
  FldArrayI posnormal(nprojectedZones,3);
  char* var; 
  E_Int m; E_Int err = 0;
  for (E_Int v  = 0 ; v < nvars; v++)
  {
    E_Int* posN = posnormal.begin(v+1);
    PyObject* l = PyList_GetItem(varsO, v);
    if (PyString_Check(l))
    {
      var = PyString_AsString(l);
      for (E_Int no = 0; no < nprojectedZones; no++)
      {
        m = K_ARRAY::isNamePresent(var, varStringP[no]);
        if (m == -1)
        {
          err = 2;
          PyErr_SetString(PyExc_TypeError,
                          "projectAllDirs: variable not found in projected zones.");
          break;
        }
        else posN[no] = m+1;
      }
    }
#if PY_VERSION_HEX >= 0x03000000
    else if (PyUnicode_Check(l))
    {
      var = (char*)PyUnicode_AsUTF8(l);
      for (E_Int no = 0; no < nprojectedZones; no++)
      {
        m = K_ARRAY::isNamePresent(var, varStringP[no]);
        if (m == -1)
        {
          err = 2;
          PyErr_SetString(PyExc_TypeError,
                          "projectAllDirs: variable not found in projected zones.");
          break;
        }
        else posN[no] = m+1;
      } 
    }
#endif 
  
    else
    {
      err = 1;
      PyErr_SetString(PyExc_TypeError,
                    "projectAllDirs: invalid string for normal component.");
    }
    
    
    if (err != 0)
    {
      for (E_Int no = 0; no < nprojectedZones; no++)
        RELEASESHAREDA(resl[no],objsP[no],fieldsP[no],a2[no],a3[no],a4[no]);   
      return NULL;
    }
  }
    
  // Extract infos from projection surfaces
  vector<E_Int> ress;  vector<char*> varStringS;
  vector<FldArrayF*> fieldsS;
  vector<void*> a2s; //ni,nj,nk ou cnt en NS
  vector<void*> a3s; //eltType en NS
  vector<void*> a4s;
  vector<PyObject*> objsS;
  skipStructured = true;
  isOk = K_ARRAY::getFromArrays(
    surfArrays, ress, varStringS, fieldsS, a2s, a3s, a4s, objsS,  
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  E_Int nsurfaces = fieldsS.size();
  if (isOk == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "projectAllDirs: invalid list of surface arrays.");
    for (E_Int no = 0; no < nprojectedZones; no++)
      RELEASESHAREDA(resl[no],objsP[no],fieldsP[no],a2[no],a3[no],a4[no]);   
    for (E_Int no = 0; no < nsurfaces; no++)
      RELEASESHAREDA(ress[no],objsS[no],fieldsS[no],a2s[no],a3s[no],a4s[no]);
    return NULL;
  }
  if (nsurfaces == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "projectAllDirs: no valid projected zone found.");
    for (E_Int no = 0; no < nprojectedZones; no++)
      RELEASESHAREDA(resl[no],objsP[no],fieldsP[no],a2[no],a3[no],a4[no]);   
    for (E_Int no = 0; no < nsurfaces; no++)
      RELEASESHAREDA(ress[no],objsS[no],fieldsS[no],a2s[no],a3s[no],a4s[no]);
    return NULL;
  }
  for (E_Int nos = 0; nos < nsurfaces; nos++)
  {
    char* eltType2 = (char*)a3s[nos];
    if (K_STRING::cmp(eltType2, "TRI") != 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "projectAllDirs: surface arrays must be TRI.");
      for (E_Int no = 0; no < nprojectedZones; no++)
        RELEASESHAREDA(resl[no],objsP[no],fieldsP[no],a2[no],a3[no],a4[no]);   
      for (E_Int no = 0; no < nsurfaces; no++)
        RELEASESHAREDA(ress[no],objsS[no],fieldsS[no],a2s[no],a3s[no],a4s[no]);
      return NULL;
    }
  }
  
  vector<E_Int> posxs; vector<E_Int> posys; vector<E_Int> poszs;
  for (E_Int no = 0; no < nsurfaces; no++)
  {
    posx1 = K_ARRAY::isCoordinateXPresent(varStringS[no]); posx1++;
    posy1 = K_ARRAY::isCoordinateYPresent(varStringS[no]); posy1++;
    posz1 = K_ARRAY::isCoordinateZPresent(varStringS[no]); posz1++;
    posxs.push_back(posx1); posys.push_back(posy1); poszs.push_back(posz1); 
  }
  
  // Build arrays
  PyObject* l = PyList_New(0);  
  vector<E_Float*> coordxp(nprojectedZones);
  vector<E_Float*> coordyp(nprojectedZones);
  vector<E_Float*> coordzp(nprojectedZones);
  PyObject* tpl;
  for (E_Int nop = 0; nop < nprojectedZones; nop++)
  {
    E_Int nfld = fieldsP[nop]->getNfld(); E_Int npts = fieldsP[nop]->getSize();
    if (resl[nop] == 1)
    {
      E_Int nip = *(E_Int*)a2[nop];
      E_Int njp = *(E_Int*)a3[nop];
      E_Int nkp = *(E_Int*)a4[nop];
      tpl = K_ARRAY::buildArray(nfld, varStringP[nop], nip, njp, nkp);
    }
    else // 2
    {
      FldArrayI* cn = (FldArrayI*)a2[nop];
      char* eltType = (char*)a3[nop];
      E_Int nelts = cn->getSize(); E_Int nvert = cn->getNfld();    
      E_Int csize = nelts*nvert;
      tpl = K_ARRAY::buildArray(nfld,varStringP[nop], npts, 
                                nelts, -1, eltType, false, csize);
      E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
      K_KCORE::memcpy__(cnnp, cn->begin(), csize);            
    }
    E_Float* fp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF f(npts, nfld, fp, true); f = *fieldsP[nop];
    coordxp[nop] = f.begin(posxp[nop]);
    coordyp[nop] = f.begin(posyp[nop]);
    coordzp[nop] = f.begin(poszp[nop]);
    PyList_Append(l, tpl); Py_DECREF(tpl);
  }

  // build BBTREE 
  typedef K_SEARCH::BoundingBox<3>  BBox3DType; 
  E_Float minB[3]; E_Float maxB[3];
  vector< vector<BBox3DType*> > vectOfBoxes(nsurfaces);// a detruire a la fin
  vector<K_SEARCH::BbTree3D*> vectOfBBTrees(nsurfaces);

  for (E_Int nos = 0; nos < nsurfaces; nos++)
  {
    FldArrayI* cnloc = (FldArrayI*)a2s[nos];
    E_Int nelts = cnloc->getSize();
    vector<BBox3DType*> boxes(nelts);// liste des bbox de ts les elements de la paroi courante
    FldArrayF bbox(nelts, 6);// xmin, ymin, zmin, xmax, ymax, zmax
    E_Float* xs = fieldsS[nos]->begin(posxs[nos]);
    E_Float* ys = fieldsS[nos]->begin(posys[nos]);
    E_Float* zs = fieldsS[nos]->begin(poszs[nos]);
    K_COMPGEOM::boundingBoxOfUnstrCells(*cnloc, xs, ys, zs, bbox);
    E_Float* xminp = bbox.begin(1); E_Float* xmaxp = bbox.begin(4);
    E_Float* yminp = bbox.begin(2); E_Float* ymaxp = bbox.begin(5);
    E_Float* zminp = bbox.begin(3); E_Float* zmaxp = bbox.begin(6);
    for (E_Int et = 0; et < nelts; et++)
    {
      minB[0] = xminp[et]; minB[1] = yminp[et]; minB[2] = zminp[et];
      maxB[0] = xmaxp[et]; maxB[1] = ymaxp[et]; maxB[2] = zmaxp[et]; 
      boxes[et] = new BBox3DType(minB, maxB);
    }
    vectOfBoxes[nos] = boxes;
    K_SEARCH::BbTree3D* bbtree = new K_SEARCH::BbTree3D(boxes);
    vectOfBBTrees[nos] = bbtree;
  }

  // Projection
  E_Float x, y, z, xo, yo, zo, dist2, distl;
  vector<E_Int> indicesBB; 
  E_Float pr1[3]; E_Float pr2[3];
  E_Float dirx, diry, dirz;
  E_Float tol = K_CONST::E_GEOM_CUTOFF;
  E_Float xsav, ysav, zsav;
  for (E_Int nop = 0; nop < nprojectedZones; nop++)
  {
    FldArrayF& fields = *fieldsP[nop];
    //E_Float* xt = fields.begin(posxp[nop]);
    //E_Float* yt = fields.begin(posyp[nop]);
    //E_Float* zt = fields.begin(poszp[nop]);
    E_Int posnx = posnormal(nop,1);
    E_Int posny = posnormal(nop,2);
    E_Int posnz = posnormal(nop,3);
    E_Float* nxt = fields.begin(posnx);
    E_Float* nyt = fields.begin(posny);
    E_Float* nzt = fields.begin(posnz);    
    E_Int npts = fields.getSize();
    E_Float* xp = coordxp[nop];
    E_Float* yp = coordyp[nop];
    E_Float* zp = coordzp[nop];

    for (E_Int ind = 0; ind < npts; ind++)
    {
      x = xp[ind]; y = yp[ind]; z = zp[ind];
      dirx = nxt[ind]; diry = nyt[ind]; dirz = nzt[ind];
      pr1[0]=x; pr2[0]=x+dirx; 
      pr1[1]=y; pr2[1]=y+diry; 
      pr1[2]=z; pr2[2]=z+dirz;
      // parcours de toutes les surfaces de projection - on conserve ensuite le projete le + proche
      dist2 = K_CONST::E_MAX_FLOAT;
      xsav = x; ysav = y; zsav = z; 
      for (E_Int nos = 0; nos < nsurfaces; nos++)
      {
        K_SEARCH::BbTree3D* bbtree = vectOfBBTrees[nos];
        bbtree->getIntersectingBoxes(pr1, pr2, indicesBB, tol);

        E_Float* xs = fieldsS[nos]->begin(posxs[nos]);
        E_Float* ys = fieldsS[nos]->begin(posys[nos]);
        E_Float* zs = fieldsS[nos]->begin(poszs[nos]);
        FldArrayI* cns = (FldArrayI*)a2s[nos];
        
        // Precond : indices : triangles candidats a la projection        
        E_Int ok = K_COMPGEOM::projectDir(x, y, z, dirx, diry, dirz,
                                          xs, ys, zs, indicesBB, *cns, xo, yo, zo, oriented);
        indicesBB.clear(); 
        if ( ok > -1 )
        {
          distl = (x-xo)*(x-xo)+(y-yo)*(y-yo)+(z-zo)*(z-zo);
          if ( distl < dist2 ) 
          {
            dist2 = distl;
            xsav = xo; ysav = yo; zsav = zo;            
          }
        }
      }
      xp[ind] = xsav; yp[ind] = ysav; zp[ind] = zsav;
    }
  }
  // Cleaning
  E_Int nboxes = vectOfBoxes.size();
  for (E_Int v0 = 0; v0 < nboxes; v0++)
  {
    vector<BBox3DType*>& boxes = vectOfBoxes[v0];
    E_Int size = boxes.size();
    for (E_Int v = 0; v < size; v++) delete boxes[v];
    delete vectOfBBTrees[v0];
  }
  vectOfBoxes.clear(); vectOfBBTrees.clear();

  for (E_Int no = 0; no < nprojectedZones; no++)
    RELEASESHAREDA(resl[no],objsP[no],fieldsP[no],a2[no],a3[no],a4[no]);   
  for (E_Int no = 0; no < nsurfaces; no++)
    RELEASESHAREDA(ress[no],objsS[no],fieldsS[no],a2s[no],a3s[no],a4s[no]);   
  return l;
}
