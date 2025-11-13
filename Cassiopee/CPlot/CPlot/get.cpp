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

#include "cplot.h"
#include "Data.h"
E_Int findNearestPoint(double xp, double yp, double zp,
                       StructZone* zone, E_Int& ind, double& dist);
E_Int findNearestPoint(double xp, double yp, double zp,
                       UnstructZone* zone, E_Int& ind, double& dist);
E_Int findElement(double xp, double yp, double zp,
                  UnstructZone* zone, double& dist, E_Int& ncon);
E_Int findElement(double xp, double yp, double zp,
                  StructZone* zone, double& dist);
E_Int findFace(double xp, double yp, double zp, E_Int elt, 
               UnstructZone* zone, double& dist);

//======================================================
// Get mode from PyObject, return an int (0: mesh, 1: solid, 
// 2: render, 3: scalar, 4: vector)
//======================================================
E_Int getMode(PyObject* modeObject)
{
  E_Int mode = -1;
  if (PyLong_Check(modeObject) || PyInt_Check(modeObject)) 
    mode = PyLong_AsLong(modeObject);
  else if (PyFloat_Check(modeObject))
    mode = int(PyFloat_AsDouble(modeObject));
  else if (PyString_Check(modeObject)) 
  {  
    char* m = PyString_AsString(modeObject);
    if (strcmp(m, "mesh")==0 || strcmp(m,"MESH")==0 || strcmp(m,"Mesh")==0) mode=0;
    if (strcmp(m, "solid")==0 || strcmp(m,"SOLID")==0 || strcmp(m,"Solid")==0) mode=1;
    if (strcmp(m, "render")==0 || strcmp(m,"RENDER")==0 || strcmp(m,"Render")==0) mode=2;
    if (strcmp(m, "scalar")==0 || strcmp(m,"SCALAR")==0 || strcmp(m,"Scalar")==0) mode=3;
    if (strcmp(m, "vector")==0 || strcmp(m,"VECTOR")==0 || strcmp(m,"Vector")==0) mode=4;
  }
#if PY_VERSION_HEX >= 0x03000000
  else if (PyUnicode_Check(modeObject))
  {
    const char* m = PyUnicode_AsUTF8(modeObject);
    if (strcmp(m, "mesh")==0 || strcmp(m,"MESH")==0 || strcmp(m,"Mesh")==0) mode=0;
    if (strcmp(m, "solid")==0 || strcmp(m,"SOLID")==0 || strcmp(m,"Solid")==0) mode=1;
    if (strcmp(m, "render")==0 || strcmp(m,"RENDER")==0 || strcmp(m,"Render")==0) mode=2;
    if (strcmp(m, "scalar")==0 || strcmp(m,"SCALAR")==0 || strcmp(m,"Scalar")==0) mode=3;
    if (strcmp(m, "vector")==0 || strcmp(m,"VECTOR")==0 || strcmp(m,"Vector")==0) mode=4;
  }
#endif  
  return mode;
}

//====================================================================
// Get scalarField from PyObject, return an int (numero de la variable)
//====================================================================
E_Int getScalarField(PyObject* scalarFieldObject)
{
  Data* d = Data::getInstance();  
  E_Int scalarField = -1;
  if (PyLong_Check(scalarFieldObject) || PyInt_Check(scalarFieldObject)) 
    scalarField = PyLong_AsLong(scalarFieldObject);
  else if (PyFloat_Check(scalarFieldObject))
    scalarField = int(PyFloat_AsDouble(scalarFieldObject));
  else if (PyString_Check(scalarFieldObject)) 
  {  
    char* m1 = PyString_AsString(scalarFieldObject);
    char m [MAXSTRINGLENGTH];
    strcpy(m, m1);
    E_Int l = strlen(m);
    if (l > 6 && m[0] == 'n' && m[1] == 'o' && m[2] == 'd' && m[3] == 'e' && 
        m[4] == 's' && m[5] == ':')
    {
      for (E_Int i = 6; i < l; i++) m[i-6] = m[i];
      m[l-6] = '\0';
    }
    if (l > 8 && m[0] == 'c' && m[1] == 'e' && m[2] == 'n' && m[3] == 't' && 
        m[4] == 'e' && m[5] == 'r' && m[6] == 's' && m[7] == ':')
    {
      for (E_Int i = 8; i < l; i++) m[i-8] = m[i];
      m[l-8] = '\0';
    }
    // cherche dans la premiere zone (if any)
    if (d->_numberOfZones > 0)
    {
      Zone* z = d->_zones[0];
      for (E_Int i = 0; i < z->nfield; i++)
      {
        if (strcmp(z->varnames[i], m) == 0) {scalarField = i; break;}
      }
    }
  }
  
#if PY_VERSION_HEX >= 0x03000000
  else if (PyUnicode_Check(scalarFieldObject))
  {
    char* m1 = (char*)PyUnicode_AsUTF8(scalarFieldObject);
    char m [MAXSTRINGLENGTH];
    strcpy(m, m1);
    E_Int l = strlen(m);
    if (l > 6 && m[0] == 'n' && m[1] == 'o' && m[2] == 'd' && m[3] == 'e' && 
        m[4] == 's' && m[5] == ':')
    {
      for (E_Int i = 6; i < l; i++) m[i-6] = m[i];
      m[l-6] = '\0';
    }
    if (l > 8 && m[0] == 'c' && m[1] == 'e' && m[2] == 'n' && m[3] == 't' && 
        m[4] == 'e' && m[5] == 'r' && m[6] == 's' && m[7] == ':')
    {
      for (E_Int i = 8; i < l; i++) m[i-8] = m[i];
      m[l-8] = '\0';
    }
    // cherche dans la premiere zone (if any)
    if (d->_numberOfZones > 0)
    {
      Zone* z = d->_zones[0];
      for (E_Int i = 0; i < z->nfield; i++)
      {
        if (strcmp(z->varnames[i], m) == 0) {scalarField = i; break;}
      }
    }
  }
#endif
  return scalarField;
}

//=============================================================================
/* 
   Get state values
*/
//=============================================================================
PyObject* K_CPLOT::getState(PyObject* self, PyObject* args)
{
  char* mode;
  if (!PyArg_ParseTuple(args, "s", &mode)) return NULL;

  Data* d = Data::getInstance();
  if (K_STRING::cmp(mode, "dim") == 0)
    return Py_BuildValue("i", d->ptrState->dim);
  else if (K_STRING::cmp(mode, "mode") == 0)
    return Py_BuildValue("i", d->ptrState->mode);
  else if (K_STRING::cmp(mode, "displayBB") == 0)
    return Py_BuildValue("i", d->ptrState->bb);
  else if (K_STRING::cmp(mode, "displayInfo") == 0)
    return Py_BuildValue("i", d->ptrState->info);
  else if (K_STRING::cmp(mode, "displayIsoLegend") == 0)
    return Py_BuildValue("i", d->ptrState->isoLegend);
  else if (K_STRING::cmp(mode, "meshStyle") == 0)
    return Py_BuildValue("i", d->ptrState->meshStyle);
  else if (K_STRING::cmp(mode, "solidStyle") == 0)
    return Py_BuildValue("i", d->ptrState->solidStyle);
  else if (K_STRING::cmp(mode, "colormap") == 0)
    return Py_BuildValue("i", d->ptrState->colormap);
  else if (K_STRING::cmp(mode, "scalarField") == 0)
    return Py_BuildValue("i", d->ptrState->scalarField);
  else if (K_STRING::cmp(mode, "scalarStyle") == 0)
    return Py_BuildValue("i", d->ptrState->scalarStyle);
  else if (K_STRING::cmp(mode, "vectorField") == 0)
    return Py_BuildValue("(i,i,i)", d->ptrState->vectorField1, 
                         d->ptrState->vectorField2, d->ptrState->vectorField3);
  else if (K_STRING::cmp(mode, "vectorStyle") == 0)
    return Py_BuildValue("i", d->ptrState->vectorStyle);
  else if (K_STRING::cmp(mode, "win") == 0)
    return Py_BuildValue("(ii)", d->_view.w, d->_view.h);
  else if (K_STRING::cmp(mode, "posCam") == 0)
    return Py_BuildValue("(fff)", d->_view.xcam,
                         d->_view.ycam, d->_view.zcam);
  else if (K_STRING::cmp(mode, "posEye") == 0)
    return Py_BuildValue("(fff)", d->_view.xeye,
                         d->_view.yeye, d->_view.zeye);
  else if (K_STRING::cmp(mode, "dirCam") == 0)
    return Py_BuildValue("(fff)", d->_view.dirx,
                         d->_view.diry, d->_view.dirz);
  else if (K_STRING::cmp(mode, "viewAngle") == 0)
    return Py_BuildValue("f", d->_view.angle);
  else if (K_STRING::cmp(mode, "ghostifyDeactivatedZones") == 0)
    return Py_BuildValue("i", d->ptrState->ghostifyDeactivatedZones);
  else if (K_STRING::cmp(mode, "niso") == 0)
    return Py_BuildValue("i", d->ptrState->niso);
  else if (K_STRING::cmp(mode, "isoEdges") == 0)
    return Py_BuildValue("f", d->ptrState->isoEdges);
  else if (K_STRING::cmp(mode, "blanking") == 0) // automatic blanking with cellN
    return Py_BuildValue("i", d->ptrState->autoblank);
  else if (K_STRING::cmp(mode, "isoScale") == 0) // iso scale of current field
  {
    E_Int nofield = d->ptrState->scalarField;
    if (nofield >= 0 && nofield < d->_nfield)
    {
      E_Int niso = d->_niso[nofield];
      double isoMin = d->_isoMin[nofield];
      double isoMax = d->_isoMax[nofield];
      if (niso == -1) 
      {
        niso = d->ptrState->niso;
        isoMin = d->minf[nofield];
        isoMax = d->maxf[nofield];
      }
      return Py_BuildValue("[idd]",niso,isoMin,isoMax);
    }
    else return Py_BuildValue("[idd]",0,0.,0.);
  }
  else if (K_STRING::cmp(mode, "colormapC1") == 0) // colormap-color1
    return Py_BuildValue("(ddd)", d->ptrState->colormapR1, d->ptrState->colormapG1, d->ptrState->colormapB1);
  else if (K_STRING::cmp(mode, "colormapC2") == 0) // colormap-color2
    return Py_BuildValue("(ddd)", d->ptrState->colormapR2, d->ptrState->colormapG2, d->ptrState->colormapB2);
  else if (K_STRING::cmp(mode, "colormapC3") == 0) // colormap-color3
    return Py_BuildValue("(ddd)", d->ptrState->colormapR3, d->ptrState->colormapG3, d->ptrState->colormapB3);
  else if (K_STRING::cmp(mode, "colormapC") == 0) // colormap-array
  {
    E_Int size = d->ptrState->colormapSize;
    float* pr = d->ptrState->colormapR;
    float* pg = d->ptrState->colormapG;
    float* pb = d->ptrState->colormapB;
    PyObject* l = PyList_New(3*size);
    for (E_Int i = 0; i < size; i++)
    {
      PyList_SetItem(l, 3*i,   Py_BuildValue("f",pr[i]));
      PyList_SetItem(l, 3*i+1, Py_BuildValue("f",pg[i]));
      PyList_SetItem(l, 3*i+2, Py_BuildValue("f",pb[i]));
    }
    return l;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "getState: unknown state field.");
    return NULL;
  }
}

//=============================================================================
/* Return selected zone 
   Retourne la derniere zone selectionnee
   Retourne -1 si aucune zone n'est actuellement selectionnee
   Sinon retourne le no de la zone selectionnee
 */
//=============================================================================
PyObject* K_CPLOT::getSelectedZone(PyObject* self, PyObject* args)
{
  Data* d = Data::getInstance();
  E_Int nz = d->ptrState->selectedZone-1;
  return Py_BuildValue("i", nz);
}

//=============================================================================
/* Return selected zones
   Retourne la liste des zones selectionnees
   Retourne une liste vide si aucune zone n'est selectionnee
*/
//=============================================================================
PyObject* K_CPLOT::getSelectedZones(PyObject* self, PyObject* args)
{
  Data* d = Data::getInstance();
  PyObject* o = PyList_New(0);
  for (E_Int i = 0; i < d->_numberOfZones; i++)
  {
    if (d->_zones[i]->selected == 1)
    {
      PyObject* l = Py_BuildValue("l", i);
      PyList_Append(o, l);
    }
  }
  return o;
}

//=============================================================================
/* Return active zones
   Retourne la liste des zones actives (affichees)
   Retourne une liste vide si aucune zone n'est affichee
*/
//=============================================================================
PyObject* K_CPLOT::getActiveZones(PyObject* self, PyObject* args)
{
  Data* d = Data::getInstance();
  PyObject* o = PyList_New(0);
  for (E_Int i = 0; i < d->_numberOfZones; i++)
  {
    if (d->_zones[i]->active == 1)
    {
      PyObject* l = Py_BuildValue("l", i);
      PyList_Append(o, l);
    }
  }
  return o;
}

//=============================================================================
/* 
   Get selected status of a zone
*/
//=============================================================================
PyObject* K_CPLOT::getSelectedStatus(PyObject* self, PyObject* args)
{
  int zone;
  if (!PYPARSETUPLE_(args, "i", &zone)) return NULL;

  Data* d = Data::getInstance();
  if (zone < 0 || zone >= d->_numberOfZones)
    return Py_BuildValue("i", -1);
  else
    return Py_BuildValue("i", d->_zones[zone]->selected);
}

//=============================================================================
/* 
   Get active status of a zone
*/
//=============================================================================
PyObject* K_CPLOT::getActiveStatus(PyObject* self, PyObject* args)
{
  int zone;
  if (!PYPARSETUPLE_(args, "i", &zone)) return NULL;

  Data* d = Data::getInstance();
  if (zone < 0 || zone >= d->_numberOfZones)
    return Py_BuildValue("i", -1);
  else
    return Py_BuildValue("i", d->_zones[zone]->active);
}

//=============================================================================
/* Return active (clicked) point coordinates 
   Retourne une liste vide si pas de point actif
   Sinon retourne la liste des coord.
 */
//=============================================================================
PyObject* K_CPLOT::getActivePoint(PyObject* self, PyObject* args)
{
  Data* d = Data::getInstance();
  E_Int nz = d->ptrState->selectedZone;
  PyObject* tpl;
  PyObject*l = PyList_New(0);
  if (nz == 0) return l;
  else
  {
    tpl = Py_BuildValue("d", d->ptrState->activePointX);
    PyList_Append(l, tpl);
    tpl = Py_BuildValue("d", d->ptrState->activePointY);
    PyList_Append(l, tpl);
    tpl = Py_BuildValue("d", d->ptrState->activePointZ);
    PyList_Append(l, tpl);
    return l;
  }
}

//=============================================================================
/* Return active (clicked) point field values 
   Retourne une liste vide si pas de point actif
   Sinon retourne la liste des valeurs des champs stockes dans le plotter.
 */
//=============================================================================
PyObject* K_CPLOT::getActivePointF(PyObject* self, PyObject* args)
{
  Data* d = Data::getInstance();
  E_Int nz = d->ptrState->selectedZone;
  PyObject* tpl;
  PyObject*l = PyList_New(0);
  if (nz == 0) return l;
  else
  {
    Zone* z = d->_zones[nz-1];
    double* activePointF = d->ptrState->activePointF;
    if (activePointF == NULL) return l;
    for (E_Int i = 0; i < z->nfield; i++)
    {
      tpl = Py_BuildValue("d", activePointF[i]);
      PyList_Append(l, tpl);
    }
    return l;
  }
}

//=============================================================================
/* Return active (clicked) point indices 
   Retourne une liste vide si pas de point actif
   Sinon retourne la liste des indices du point actif.
 */
//=============================================================================
PyObject* K_CPLOT::getActivePointIndex(PyObject* self, PyObject* args)
{
  Data* d = Data::getInstance();
  E_Int nz = d->ptrState->selectedZone;
  PyObject* tpl;
  PyObject*l = PyList_New(0);
  if (nz == 0) return l;
  else
  {
    double posX, posY, posZ;
    E_Int zone, ind, indE, ncon; double dist;
    posX = d->ptrState->activePointX;
    posY = d->ptrState->activePointY;
    posZ = d->ptrState->activePointZ;

    // recompute zone and index
    //d->findBlockContaining(posX, posY, posZ, 
    //                       zone, ind, indE, dist, ncon);

    // only recompute index
    double dn, de;
    zone = d->ptrState->selectedZone-1;
    Zone* z = d->_zones[zone];
    if (zone < d->_numberOfStructZones)
    {
      findNearestPoint(posX, posY, posZ, (StructZone*)z, ind, dn);
      indE = findElement(posX, posY, posZ, (StructZone*)z, de);
    }
    else
    {
      findNearestPoint(posX, posY, posZ, (UnstructZone*)z, ind, dn);
      indE = findElement(posX, posY, posZ, (UnstructZone*)z, de, ncon);
    }
    
    if (zone < d->_numberOfStructZones)
    {
      StructZone* zz = (StructZone*)z;
      E_Int ni = zz->ni; 
      E_Int nj = zz->nj;
      E_Int k = ind / (ni*nj);
      E_Int j = (ind - k*ni*nj)/ni;
      E_Int i = ind - k*ni*nj - j*ni;
      d->ptrState->activePointI = i+1;
      d->ptrState->activePointJ = j+1;
      d->ptrState->activePointK = k+1;
    }
    else
    {
      d->ptrState->activePointI = ind; // indice du noeud le plus proche
      d->ptrState->activePointJ = indE; // indice de l'element contenant P
      d->ptrState->activePointL = ncon; // connectivite contenant l'element
      UnstructZone* zz = (UnstructZone*)z;
      if (zz->eltType[0] != 10) // autre que NGON
      {
        E_Int* c = zz->connect[ncon];
        E_Int size = zz->eltSize[ncon];
        E_Int ne = zz->nec[ncon];
        E_Int v = 0;
        E_Int prev = 0;
        for (E_Int nc = 0; nc < ncon; nc++) prev += zz->nec[nc];
        for (v = 0; v < size; v++)
        {
          if (c[indE-prev+v*ne] == ind+1) break;
        }
        d->ptrState->activePointK = -v-1;
      }
      else d->ptrState->activePointK = -findFace(
        posX, posY, posZ, indE, zz, dist);
    }

    if (d->ptrState->activePointK < 0) // non structure
    {
      // indice noeud le plus proche
      tpl = Py_BuildValue("l", d->ptrState->activePointI);
      PyList_Append(l, tpl);
      // indice centre le plus proche
      tpl = Py_BuildValue("l", d->ptrState->activePointJ);
      PyList_Append(l, tpl);
      // 3 zeros 
      tpl = Py_BuildValue("l", abs(d->ptrState->activePointK));
      PyList_Append(l, tpl);
      tpl = Py_BuildValue("l", 0);
      PyList_Append(l, tpl);
      tpl = Py_BuildValue("l", 0);
      PyList_Append(l, tpl);
    }
    else
    { // structure
      StructZone* z = (StructZone*)d->_zones[nz-1];
      E_Int ni = z->ni;
      E_Int nj = z->nj;
      E_Int indi = d->ptrState->activePointI;
      E_Int indj = d->ptrState->activePointJ;
      E_Int indk = d->ptrState->activePointK;
      E_Int ind = indi-1 + (indj-1)*ni + (indk-1)*ni*nj;
      tpl = Py_BuildValue("l", ind);
      PyList_Append(l, tpl);
      E_Int indc = indi-1 + (indj-1)*(ni-1) + (indk-1)*(ni-1)*(nj-1);
      tpl = Py_BuildValue("l", indc);
      PyList_Append(l, tpl);
      tpl = Py_BuildValue("l", indi);
      PyList_Append(l, tpl);
      tpl = Py_BuildValue("l", indj);
      PyList_Append(l, tpl);
      tpl = Py_BuildValue("l", indk);
      PyList_Append(l, tpl);
    }
    return l;
  }
}

//=============================================================================
/* 
   Return keyboard. 
   Retourne la chaine du keyboard.
 */
//=============================================================================
PyObject* K_CPLOT::getKeyboard(PyObject* self, PyObject* args)
{
  Data* d = Data::getInstance();
  PyObject* tpl;
  char tmp[128];
  for (int i = 0; i < d->ptrState->kcursor; i++)
  {
    tmp[i] = d->ptrState->keys[i];
  }
  tmp[d->ptrState->kcursor] = '\0';
  tpl = Py_BuildValue("s", tmp);
  return tpl;

  // modification pour forcer le latin1
  /*
  int j = 0;
  for (int i = 0; i < d->ptrState->kcursor; i++)
  {
    unsigned char c = d->ptrState->keys[i];
    unsigned int k = (unsigned int)(c);
    printf("%u %u\n", c, k);
    //if (k == -61) { tmp[j] = 0x00; tmp[j+1] = 0xe9; j+=2; }
    if (k == -61) { tmp[j] = 233; j+=1; }
    else { tmp[j] = c; j+=1; }
  }
  tmp[j] = '\0';
  tpl = PyUnicode_DecodeLatin1(tmp, strlen(tmp), "ignore");
  return tpl;
  */
}

//=============================================================================
/*
  Return mouse position (coordinates in 3D space) and the mouse button
  state while dragging (GLUT_MIDDLE_BUTTON(1), GLUT_LEFT_BUTTON (0), 
  GLUT_RIGHT_BUTTON (2), RELEASED (5))
*/
//=============================================================================
PyObject* K_CPLOT::getMouseState(PyObject* self, PyObject* args)
{
  Data* d = Data::getInstance();
  float posX = d->ptrState->currentMousePosX;
  float posY = d->ptrState->currentMousePosY;
  float posZ = d->ptrState->currentMousePosZ;
  E_Int mouseButton = d->ptrState->currentMouseButton;
  //printf("%f %f %f\n", posX, posY, posZ);
  PyObject* tpl = Py_BuildValue("lddd", mouseButton, posX, posY, posZ);
  return tpl;
}
