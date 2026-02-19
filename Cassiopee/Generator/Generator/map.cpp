/*    
    Copyright 2013-2026 ONERA.

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

// map distribution on a curve

# include "generator.h"
# include <vector>
# include "CompGeom/compGeom.h"
using namespace std; 
using namespace K_FLD;

//=============================================================================
// Map a 1D-distribution on a curve
//=============================================================================
PyObject* K_GENERATOR::mapMesh( PyObject* self, PyObject* args )
{
  PyObject* array; PyObject* arrayd;
  if (!PYPARSETUPLE_(args, OO_, &array, &arrayd)) return NULL;

  // Check array
  E_Int ni, nj, nk, nid, njd, nkd;
  FldArrayF* f; FldArrayF* fd;
  char* varString; char* eltType;
  FldArrayI* cn; FldArrayI* cnd;
  char* varStringd; char* eltTyped;

  // Extraction des infos sur le maillage
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, cn, eltType);
  // Extraction des infos sur la distribution
  E_Int resd = K_ARRAY::getFromArray3(arrayd, varStringd, fd,
                                      nid, njd, nkd, cnd, eltTyped);

  E_Int api = f->getApi();
  if (res == 1 && resd == 1)
  {
    // find x,y,z if possible
    E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
    E_Int posxd = K_ARRAY::isCoordinateXPresent(varStringd);
    E_Int posyd = K_ARRAY::isCoordinateYPresent(varStringd);
    E_Int poszd = K_ARRAY::isCoordinateZPresent(varStringd);
  
    if (posx == -1 || posy == -1 || posz == -1 ||
        posxd == -1 || posyd == -1 || poszd == -1)
    {
      RELEASESHAREDB(res, array, f, cn);
      RELEASESHAREDB(resd, arrayd, fd, cnd);
      PyErr_SetString(PyExc_TypeError,
                      "map: coordinates not found.");
      return NULL;
    }
    posx++; posy++; posz++;
    posxd++; posyd++; poszd++;
    
    // Mapping sur courbe
    FldArrayF s(ni);
    FldArrayF dx(ni);
    FldArrayF dy(ni);
    FldArrayF dz(ni);
    
    PyObject* tpl = K_ARRAY::buildArray3(3, "x,y,z", nid, 1, 1, api);
    FldArrayF* coord1;
    K_ARRAY::getFromArray3(tpl, coord1);

    K_COMPGEOM::onedmap(
      ni, f->begin(posx), f->begin(posy), f->begin(posz),
      nid, fd->begin(posxd),
      coord1->begin(1), coord1->begin(2), coord1->begin(3),
      s.begin(), dx.begin(), dy.begin(), dz.begin()
    );

    s.malloc(0); dx.malloc(0); dy.malloc(0); dz.malloc(0);
    RELEASESHAREDS(tpl, coord1);
    RELEASESHAREDB(res, array, f, cn);
    RELEASESHAREDB(resd, arrayd, fd, cnd);
    return tpl;
  }
  else if (res == 2 && resd == 1)
  {
    // find x,y,z if possible
    E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
    E_Int posxd = K_ARRAY::isCoordinateXPresent(varStringd);
    E_Int posyd = K_ARRAY::isCoordinateYPresent(varStringd);
    E_Int poszd = K_ARRAY::isCoordinateZPresent(varStringd);
  
    if (posx == -1 || posy == -1 || posz == -1 ||
        posxd == -1 || posyd == -1 || poszd == -1)
    {
      RELEASESHAREDB(res, array, f, cn);
      RELEASESHAREDB(resd, arrayd, fd, cnd);
      PyErr_SetString(PyExc_TypeError,
                      "map: coordinates not found.");
      return NULL;
    }
    posx++; posy++; posz++;
    posxd++; posyd++; poszd++;
    E_Int ni = f->getSize();
    E_Int net0 = cn->getSize();
    E_Int* cn1 = cn->begin(1);
    E_Int* cn2 = cn->begin(2);
    // Mapping sur courbe
    FldArrayF s(ni); s.setAllValuesAtNull(); // tableaux de travail
    FldArrayF dx(ni); dx.setAllValuesAtNull();
    FldArrayF dy(ni); dy.setAllValuesAtNull();
    FldArrayF dz(ni); dz.setAllValuesAtNull();
   
    // connectivite sortante: verif bar entrante fermee
    E_Int net = nid-1;
    E_Int inds = cn1[0]; E_Int inde = cn2[net0-1];
    //if (inds == inde) net = net-1;

    PyObject* tpl = K_ARRAY::buildArray3(3, "x,y,z", nid, net,
                                         "BAR", false, api);
    FldArrayF* coord1; FldArrayI* cnout;
    K_ARRAY::getFromArray3(tpl, coord1, cnout);

    K_COMPGEOM::onedmapbar(
      ni, f->begin(posx), f->begin(posy), f->begin(posz),
      nid, fd->begin(posxd), net0, cn1, cn2, 
      net, cnout->begin(1), cnout->begin(2),
      coord1->begin(1), coord1->begin(2), coord1->begin(3),
      s.begin(), dx.begin(), dy.begin(), dz.begin()
    );
    s.malloc(0); dx.malloc(0); dy.malloc(0); dz.malloc(0);

    if (inds == inde) { (*cnout)(net-1, 2) = 1; }

    RELEASESHAREDU(tpl, coord1, cnout);
    RELEASESHAREDB(res, array, f, cn);
    RELEASESHAREDB(resd, arrayd, fd, cnd);
    return tpl;
  }
  else 
  {
    RELEASESHAREDB(res, array, f, cn);
    RELEASESHAREDB(resd, arrayd, fd, cnd);
    PyErr_SetString(PyExc_TypeError,
                    "map: unknown type of arrays.");
    return NULL;
  }
}
