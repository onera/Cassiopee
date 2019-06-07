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
// File conversion

#include <stdio.h>
#include <string.h>
#include "converter.h"
#include "kcore.h"
#include "IO/GenIO.h"

using namespace std;
using namespace K_FLD;

// ============================================================================
/* Convert file to arrays */
// ============================================================================
PyObject* K_CONVERTER::convertFile2Arrays(PyObject* self, PyObject* args)
{
  E_Int n;
  E_Int nLine; E_Int nCurve; E_Float density;
  char* fileName; char* fileFmt;
  PyObject* zoneNamesO; PyObject* BCFacesO;

  if (!PYPARSETUPLE(args, "sslldOO", "ssiidOO",
                    "ssllfOO", "ssiifOO", &fileName, &fileFmt, 
                    &nCurve, &nLine, &density, &zoneNamesO,
                    &BCFacesO)) return NULL;

  E_Int NptsLine = nLine;
  E_Int NptsCurve = nCurve;
  // Check recognised formats
  if (checkRecognisedFormat(fileFmt) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "convertFile2Arrays: unknown file format.");
    return NULL;
  }

  // Check zoneNames list
  E_Int znOsize;
  E_Bool buildZoneNames = true;
  if (zoneNamesO == Py_None)
  {
    buildZoneNames = false; // on ne construit la liste des noms de zone 
    znOsize = 0;
  }
  else if (PyList_Check(zoneNamesO) != 0)
    znOsize = PyList_Size(zoneNamesO);
  else
  {
    PyErr_SetString(
      PyExc_TypeError, 
      "convertFile2Arrays: zoneNames argument must be either a list or None.");
    return NULL;
  }
  if (znOsize > 0)
  {
    printf("Warning: convertFile2Arrays: zoneNames is not empty. zones names will be appended to zoneNames.\n");
  }

  vector<FldArrayF*> field;            // field read for each zone
  vector<E_Int> im; vector<E_Int> jm; vector<E_Int> km;
  char* varString = NULL;
  vector<FldArrayI*> c;
  vector<FldArrayF*> ufield;
  vector<E_Int> et;
  vector<char*> zoneNames;
  vector<FldArrayI*> BCFaces;
  vector<char*> BCNames;
  E_Int ret = 1;

  printf("Reading %s (%s)...", fileName, fileFmt);
  fflush(stdout);

  if (K_STRING::cmp(fileFmt, "bin_tp") == 0)
  {
    // Binary tecplot read
    ret = K_IO::GenIO::getInstance()->tecread(fileName, varString, field, 
                                              im, jm, km, 
                                              ufield, c, et, zoneNames);
  }
  else if (K_STRING::cmp(fileFmt, "fmt_tp") == 0)
  {
    // Formatted tecplot read
    ret = K_IO::GenIO::getInstance()->tpread(fileName, varString, field,
                                             im, jm, km, 
                                             ufield, c, et, zoneNames);
  }
  else if (K_STRING::cmp(fileFmt, "fmt_v3d") == 0)
  {
    // Formatted v3d read
    ret = K_IO::GenIO::getInstance()->fv3dread(fileName, varString, 
                                               im, jm, km, field, zoneNames);
  }
  else if (K_STRING::cmp(fileFmt, "bin_v3d") == 0)
  {
    // Binary v3d read
    ret = K_IO::GenIO::getInstance()->v3dread(fileName, varString, 
                                              im, jm, km, field, zoneNames);
  }
  else if (K_STRING::cmp(fileFmt, "fmt_plot3d") == 0)
  {
    // Formatted plot3d read
    ret = K_IO::GenIO::getInstance()->fp3dread(fileName, varString, 
                                               im, jm, km, field, zoneNames);
  }
  else if (K_STRING::cmp(fileFmt, "bin_plot3d") == 0)
  {
    // Binary plot3d
    ret = K_IO::GenIO::getInstance()->plot3dread(fileName, varString, 
                                                 im, jm, km, field, zoneNames);
  } 
  else if (K_STRING::cmp(fileFmt, "fmt_pov") == 0)
  {
    // Formatted pov read
    ret = K_IO::GenIO::getInstance()->povread(fileName, varString, field, 
                                              im, jm, km, 
                                              ufield, c, et, zoneNames);
  }
  else if (K_STRING::cmp(fileFmt, "fmt_mesh") == 0)
  {
    // Formatted mesh read
    ret = K_IO::GenIO::getInstance()->meshread(fileName, varString, field, 
                                               im, jm, km, 
                                               ufield, c, et, zoneNames);
  }
  else if (K_STRING::cmp(fileFmt, "fmt_gmsh") == 0)
  {
    // Formatted gmsh read
    ret = K_IO::GenIO::getInstance()->gmshread(fileName, varString, field, 
                                               im, jm, km, 
                                               ufield, c, et, zoneNames);
  }
  else if (K_STRING::cmp(fileFmt, "bin_gmsh") == 0)
  {
    // Formatted gmsh read
    ret = K_IO::GenIO::getInstance()->bingmshread(fileName, varString, field, 
                                                  im, jm, km, 
                                                  ufield, c, et, zoneNames);
  }
  else if (K_STRING::cmp(fileFmt, "fmt_obj") == 0)
  {
    // Formatted obj read
    ret = K_IO::GenIO::getInstance()->objread(fileName, varString, field, 
                                              im, jm, km, 
                                              ufield, c, et, zoneNames);
  }
  else if (K_STRING::cmp(fileFmt, "bin_stl") == 0)
  {
    // Binary STL read
    ret = K_IO::GenIO::getInstance()->stlread(fileName, varString, field, 
                                              im, jm, km, 
                                              ufield, c, et, zoneNames);
  }
  else if (K_STRING::cmp(fileFmt, "fmt_stl") == 0)
  {
    // Formatted STL read
    ret = K_IO::GenIO::getInstance()->fstlread(fileName, varString, field, 
                                               im, jm, km, 
                                               ufield, c, et, zoneNames);
  }
  else if (K_STRING::cmp(fileFmt, "bin_3ds") == 0)
  {
    // Binary 3DS read
    ret = K_IO::GenIO::getInstance()->f3dsread(fileName, varString, field, 
                                               im, jm, km, 
                                               ufield, c, et, zoneNames);
  }
  else if (K_STRING::cmp(fileFmt, "bin_ply") == 0)
  {
    // Binary PLY read
    ret = K_IO::GenIO::getInstance()->plyread(fileName, varString, field, 
                                              im, jm, km, 
                                              ufield, c, et, zoneNames);
  }
  else if (K_STRING::cmp(fileFmt, "fmt_xfig") == 0)
  {
    // Formatted Xfig read
    ret = K_IO::GenIO::getInstance()->xfigread(fileName, varString, 
                                               NptsCurve, NptsLine, 
                                               field, im, jm, km, 
                                               ufield, c, et, zoneNames);
  }
  else if (K_STRING::cmp(fileFmt, "fmt_svg") == 0)
  {
    // Formatted SVG read
    /*
    if (density == -1.)
      ret = K_IO::GenIO::getInstance()->svgread(fileName, varString, 
                                                NptsCurve, NptsLine, 
                                                field, im, jm, km, 
                                                ufield, c, et, zoneNames);
    else 
      ret = K_IO::GenIO::getInstance()->svgread(fileName, varString, 
                                                density,
                                                field, im, jm, km, 
                                                ufield, c, et, zoneNames);
                                                */
    ret = K_IO::GenIO::getInstance()->svgread(fileName, varString, density,
                                              NptsCurve, NptsLine, 
                                              field, im, jm, km, 
                                              ufield, c, et, zoneNames);
    
  }
  else if (K_STRING::cmp(fileFmt, "fmt_gts") == 0)
  {
    // Formatted GTS read
    ret = K_IO::GenIO::getInstance()->gtsread(fileName, varString, field, 
                                              im, jm, km, 
                                              ufield, c, et, zoneNames);
  }
  else if (K_STRING::cmp(fileFmt, "fmt_iges") == 0)
  {
    // Formatted IGES read
    ret = K_IO::GenIO::getInstance()->igesread(fileName, varString, field, 
                                               im, jm, km, 
                                               ufield, c, et, zoneNames);
  }
  else if (K_STRING::cmp(fileFmt, "bin_png") == 0)
  {
    // Binary PNG read
    ret = K_IO::GenIO::getInstance()->pngread(fileName, varString, field, 
                                              im, jm, km, 
                                              ufield, c, et, zoneNames);
  }
  else if (K_STRING::cmp(fileFmt, "fmt_su2") == 0)
  {
    // Formatted su2 read
    ret = K_IO::GenIO::getInstance()->su2read(fileName, varString, 
                                              field, im, jm, km, 
                                              ufield, c, et, zoneNames,
                                              BCFaces, BCNames);
    // Pour l'instant, on ne peut pas les traiter en elements
    for (size_t i = 0; i < BCFaces.size(); i++) delete BCFaces[i];
    for (size_t i = 0; i < BCNames.size(); i++) delete [] BCNames[i];
    BCFaces.clear(); BCNames.clear();
  }
  else if (K_STRING::cmp(fileFmt, "fmt_cedre") == 0)
  {
    // Formatted cedre read
    ret = K_IO::GenIO::getInstance()->cedreread(fileName, varString, field, 
                                                im, jm, km, 
                                                ufield, c, et, zoneNames,
                                                BCFaces, BCNames);
  }
  else if (K_STRING::cmp(fileFmt, "fmt_tgf") == 0)
  {
    // Formatted tgf read
    ret = K_IO::GenIO::getInstance()->tgfread(fileName, varString, field, 
                                              im, jm, km, 
                                              ufield, c, et, zoneNames,
                                              BCFaces, BCNames);
  }
  else
  {
    // Default: error
    PyErr_SetString(PyExc_TypeError,
                    "convertFile2Arrays: unrecognised format.");
    return NULL;
  }

  // Ajout des BCFaces dans l'arbre
  E_Int size = BCFaces.size();
  if (BCFacesO != Py_None)
  {
    char n[128]; E_Int c = 0;
    for (E_Int j = 0; j < size; j++) // pour chaque bloc
    {
      char* myBCNames = BCNames[j];
      FldArrayI& myBCFaces = *BCFaces[j];
      vector<char*> names; // unique BC names
      E_Int np = myBCFaces.getSize();
      FldArrayI indir(np);
      E_Int l, k, nbc;
      E_Boolean exist;
      
      for (E_Int i = 0; i < np; i++) // pour chaque face frontiere
      {
        // Recopie du nom de BC
        l = 0;
        while (myBCNames[c] != '\0') { n[l] = myBCNames[c]; c++; l++; }
        n[l] = '\0'; c++;
        //printf("%s\n", myBCNames);
        // BC deja existante?
        exist = false;
        for (unsigned int k = 0; k < names.size(); k++)
        {
          if (K_STRING::cmp(n, names[k]) == 0) { indir[i] = k; exist = true; break; }
        }
        if (exist == false)
        { char* na = new char[128]; strcpy(na, n); 
          names.push_back(na); indir[i] = names.size()-1; }
      }
      //for (E_Int i = 0; i < names.size(); i++) printf("%s\n", names[i]);

      // Construction de l'objet python de sortie
      PyObject* lc = PyList_New(0);
      E_Int sizeNames = names.size();
      for (E_Int i = 0; i < sizeNames; i++)
      {
        // nom de la BC
        PyObject* charObject = Py_BuildValue("s", names[i]);
        nbc = 0;
        for (E_Int j = 0; j < np; j++) { if (indir[j] == i) nbc++; }
        PyObject* a = K_NUMPY::buildNumpyArray(nbc, 1, 1);
        E_Int* pt = K_NUMPY::getNumpyPtrI(a);
        
        k = 0;
        for (E_Int j = 0; j < np; j++)
        { if (indir[j] == i) { pt[k] = myBCFaces[j]; k++; } } 
        // delete chars
        delete [] names[i];
        PyList_Append(lc, charObject); Py_DECREF(charObject);
        PyList_Append(lc, a); Py_DECREF(a);
      }
      PyList_Append(BCFacesO, lc); Py_DECREF(lc);
    }
  }
  for (E_Int i = 0; i < size; i++) delete BCFaces[i];
  for (E_Int i = 0; i < size; i++) delete [] BCNames[i];
  /* Fin BCFaces */

  if (ret == 1)
  {
    char error[256];
    sprintf(error, "convertFile2Arrays: fail to read %s (%s).", 
            fileName, fileFmt);
    PyErr_SetString(PyExc_IOError, error);
    return NULL;
  }
  printf("done.\n"); fflush(stdout);

  // Building numpy arrays
  PyObject* tpl;
    
  if (im.size() == 0 && c.size() == 0)  
  {
    printf("Warning: convertFile2Arrays: no block in file.\n");
    return PyList_New(0);
  }
  
  PyObject* l = PyList_New(0);
  
  n = im.size();
  for (E_Int i = 0; i < n; i++)
  {
    // Build array
    tpl = K_ARRAY::buildArray2(*field[i], varString,
                               im[i], jm[i], km[i]);
    delete field[i];
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }
  
  n = ufield.size();    
  for (E_Int i = 0; i < n; i++)
  {
    char eltType[28]; E_Int d;
    K_ARRAY::typeId2eltString(et[i], 0, eltType, d);
    tpl = K_ARRAY::buildArray2(*ufield[i], varString,
                               *c[i], eltType);
    delete ufield[i]; delete c[i];
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }

  // build zoneNames list. Les fonctions de lecture ont alloue un char* par
  // zone lue dans l'ordre (zones structurees, non structurees)
  if (buildZoneNames)
  {
    E_Int znsize = zoneNames.size();
    for (E_Int i = 0; i < znsize; i++)
    {
      PyObject* charObject; 
      charObject = Py_BuildValue("s", zoneNames[i]);
      PyList_Append(zoneNamesO, charObject);
      Py_DECREF(charObject);
      delete [] zoneNames[i];
    }
  }
  else
  {
    E_Int znsize = zoneNames.size();
    for (E_Int i = 0; i < znsize; i++) delete [] zoneNames[i];
  }

  delete [] varString;

  return l;
}

//============================================================================ 
/* Convert arrays to file */
//============================================================================
PyObject* K_CONVERTER::convertArrays2File(PyObject* self, PyObject* args)
{
  PyObject* arrays;
  E_Int is, rs, c;
  char* e; char* fileName; char* fileFmt; char* dataFmt;
  char* eltType;
  PyObject* zoneNamesO; PyObject* BCFacesO;

  if (!PYPARSETUPLEI(args, "OssllslsOO", "OssiisisOO", 
                     &arrays, &fileName, &fileFmt, &is, &rs, &e, &c, 
                     &dataFmt, &zoneNamesO, &BCFacesO)) return NULL;

  E_Int isize     = is;
  E_Int rsize     = rs;
  E_Int endianess = 1;
  E_Int colormap  = c;

  if (K_STRING::cmp(e, "little") == 0) endianess = 0;
  else if (K_STRING::cmp(e, "big") == 0) endianess = 1;
  else
    printf("Warning: convertArrays2File: endian must be 'little' or 'big'. Set to 'big'.\n");

  // Check recognised formats
  if (checkRecognisedFormat(fileFmt) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "convertArrays2File: unknown file format.");
    return NULL;
  }

  // Check every array
  if (PyList_Check(arrays) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "convertArrays2File: arrays argument must be a list.");
    return NULL;
  }

  // Check zonenames list
  if (PyList_Check(zoneNamesO) == 0)
  {
    PyErr_SetString(
      PyExc_TypeError, 
      "convertArrays2File: zoneNames argument must be a list.");
    return NULL;
  }
  // Building zoneNames vector
  vector<char*> zoneNames;
  for (int z = 0; z < PyList_Size(zoneNamesO); z++)
  {
    PyObject* tplz = PyList_GetItem(zoneNamesO, z);
    if (PyString_Check(tplz))
    {
      char* str = PyString_AsString(tplz);
      zoneNames.push_back(str);
    }
#if PY_VERSION_HEX >= 0x03000000
    else if (PyUnicode_Check(tplz))
    {
      char* str = PyBytes_AsString(PyUnicode_AsUTF8String(tplz));
      zoneNames.push_back(str);  
    }
#endif  
    else
    {
      printf("Warning: convertArrays2File: zone name must be a string. Skipped...\n");
    }
  }

  // Detection endian
  E_Int endianChange;
  E_Int machineEndianess = K_IO::GenIO::getInstance()->machineEndianess();  
  if (machineEndianess == endianess) endianChange = 0;
  else endianChange = 1;

  E_Int n = PyList_Size(arrays);
  
  PyObject* tpl;
  E_Int nil, njl, nkl;
  vector<E_Int> ni; vector<E_Int> nj; vector<E_Int> nk;
  vector<FldArrayF*> fieldc;
  char* varString;
  vector<E_Int> elt;
  FldArrayF* f; FldArrayI* cn;
  E_Int res;
  vector<FldArrayF*> fieldu; // non structuré
  vector<FldArrayI*> connectu;
  for (int i = 0; i < n; i++)
  {
    tpl = PyList_GetItem(arrays, i);
    res = K_ARRAY::getFromArray(tpl, varString, 
                                f, nil, njl, nkl, cn, eltType);
      
    if (res == 1)
    {
      if (nil*njl*nkl > 0)
      {
        ni.push_back(nil); nj.push_back(njl); nk.push_back(nkl);
        fieldc.push_back(f);
      }
      else 
        printf("Warning: convertArrays2File: one array is empty.\n");
    }
    else if (res == 2)
    {
      if (f->getSize() > 0)
      {
        E_Int id = getElementTypeId(eltType);
        // Ecriture non structuree
        if (id == 0) // NODE-> structure
        {
          ni.push_back(f->getSize()); nj.push_back(1); nk.push_back(1);
          fieldc.push_back(f); 
        }
        else if (id > 0)
        {
          fieldu.push_back(f);
          connectu.push_back(cn);
          elt.push_back(id);
        }
        else //id == -1
          printf("Warning: convertArrays2File: element type %s not taken into account.\n", eltType);
      }
      else 
        printf("Warning: convertArrays2File: one array is empty.\n");
    }
    else
      printf("Warning: convertArrays2File: one array is invalid.\n");
  }

  // Nfld
  if (fieldc.size() == 0 && fieldu.size() == 0)
  {
    printf("Warning: convertArrays2File: nothing to write.\n");
    Py_INCREF(Py_None);
    return Py_None;
  }
 
  // Writing output file
  printf("Writing %s (%s)...", fileName, fileFmt);
  fflush(stdout);
  E_Int isok;
  
  if (K_STRING::cmp(fileFmt, "bin_tp") == 0) // binary tecplot
  {
    isok = 
      K_IO::GenIO::getInstance()->tecwrite(fileName, dataFmt, varString,
                                           ni, nj, nk,
                                           fieldc, fieldu, connectu, elt,
                                           zoneNames);
  }
  else if (K_STRING::cmp(fileFmt, "fmt_tp") == 0) // fmt tecplot
  { 
    isok = K_IO::GenIO::getInstance()->tpwrite(fileName, dataFmt, varString,
                                               ni, nj, nk,
                                               fieldc, fieldu, connectu, elt,
                                               zoneNames);
  }
  else if (K_STRING::cmp(fileFmt, "fmt_v3d") == 0) // fmt v3d
  {
    if (fieldu.size() != 0)
      printf("Warning: convertArrays2File: unstructured arrays not converted.\n"); 
    
    isok = K_IO::GenIO::getInstance()->fv3dwrite(fileName, dataFmt, varString,
                                                 ni, nj, nk,
                                                 fieldc, zoneNames);
  }
  else if (K_STRING::cmp(fileFmt, "bin_v3d") == 0) // binary v3d
  {
    if (fieldu.size() != 0)
      printf("Warning: convertArrays2File: unstructured arrays not converted.\n"); 
    
    if (rsize == 4)
    {
      printf("Warning: convertArrays2File: r4 option is not supported.\n");
      rsize = 8;
    }
    isok = K_IO::GenIO::getInstance()->v3dwrite(fileName, dataFmt, varString,
                                                ni, nj, nk,
                                                fieldc, zoneNames, 
                                                isize, rsize, endianChange);
  }
  else if (K_STRING::cmp(fileFmt, "fmt_plot3d") == 0) // fmt plot3d
  {
    if (fieldu.size() != 0)
      printf("Warning: convertArrays2File: unstructured arrays not converted.\n"); 
    
    isok = K_IO::GenIO::getInstance()->fp3dwrite(fileName, dataFmt, varString,
                                                 ni, nj, nk,
                                                 fieldc, zoneNames);
  }
  else if (K_STRING::cmp(fileFmt, "bin_plot3d") == 0) // binary plot3d
  {
    if (fieldu.size() != 0)
      printf("Warning: convertArrays2File: unstructured arrays not converted.\n"); 
    
    isok = K_IO::GenIO::getInstance()->plot3dwrite(
      fileName, dataFmt, varString, ni, nj, nk,
      fieldc, zoneNames, isize, rsize, endianChange);
  }
  else if (K_STRING::cmp(fileFmt, "fmt_pov") == 0) // fmt pov
  {
    if (fieldc.size() != 0)
      printf("Warning: convertArrays2File: structured arrays not converted.\n"); 
    
    isok = K_IO::GenIO::getInstance()->povwrite(fileName, dataFmt, varString,
                                                ni, nj, nk,
                                                fieldc, fieldu, connectu, elt,
                                                zoneNames, colormap);
  }
  else if (K_STRING::cmp(fileFmt, "bin_df3") == 0) // binary df3
  {
    if (fieldu.size() != 0)
      printf("Warning: convertArrays2File: unstructured arrays not converted.\n"); 
    
    isok = K_IO::GenIO::getInstance()->df3write(fileName, dataFmt, varString,
                                                ni, nj, nk,
                                                fieldc, fieldu, connectu, elt,
                                                zoneNames);
  } 
  else if (K_STRING::cmp(fileFmt, "fmt_mesh") == 0) // fmt mesh
  {
    if (fieldc.size() != 0)
      printf("Warning: convertArrays2File: structured arrays not converted.\n"); 
    
    isok = K_IO::GenIO::getInstance()->meshwrite(fileName, dataFmt, varString,
                                                 ni, nj, nk,
                                                 fieldc, fieldu, connectu, elt,
                                                 zoneNames);
  }
  else if (K_STRING::cmp(fileFmt, "fmt_gmsh") == 0) // fmt gmsh
  {
    isok = K_IO::GenIO::getInstance()->gmshwrite(fileName, dataFmt, varString,
                                                 ni, nj, nk,
                                                 fieldc, fieldu, connectu, elt,
                                                 zoneNames);
  }
  else if (K_STRING::cmp(fileFmt, "bin_gmsh") == 0) // bin gmsh
  {
    isok = K_IO::GenIO::getInstance()->bingmshwrite(
      fileName, dataFmt, varString,
      ni, nj, nk,
      fieldc, fieldu, connectu, elt,
      zoneNames);
  }
  else if (K_STRING::cmp(fileFmt, "bin_png") == 0) // in png
  {
    if (fieldu.size() != 0)
      printf("Warning: convertArrays2File: unstructured arrays not converted.\n"); 
    
    isok = K_IO::GenIO::getInstance()->pngwrite(fileName, dataFmt, varString,
                                                ni, nj, nk,
                                                fieldc, fieldu, connectu, elt,
                                                zoneNames);
  }
  else if (K_STRING::cmp(fileFmt, "fmt_su2") == 0) // fmt su2
  {
    if (fieldc.size() != 0)
      printf("Warning: convertArrays2File: structured arrays not converted.\n"); 
    
    isok = K_IO::GenIO::getInstance()->su2write(fileName, dataFmt, varString,
                                                ni, nj, nk,
                                                fieldc, fieldu, connectu, elt,
                                                zoneNames, BCFacesO);
  }
  else if (K_STRING::cmp(fileFmt, "fmt_obj") == 0) // fmt obj
  {
    if (fieldc.size() != 0)
      printf("Warning: convertArrays2File: structured arrays not converted.\n"); 
    
    isok = K_IO::GenIO::getInstance()->objwrite(fileName, dataFmt, varString,
                                                ni, nj, nk,
                                                fieldc, fieldu, connectu, elt, 
                                                zoneNames);
  }
  else if (K_STRING::cmp(fileFmt, "bin_stl") == 0) // bin stl
  {
    if (fieldc.size() != 0)
      printf("Warning: convertArrays2File: structured arrays not converted.\n"); 
    
    isok = K_IO::GenIO::getInstance()->stlwrite(fileName, dataFmt, varString,
                                                ni, nj, nk,
                                                fieldc, fieldu, connectu, elt,
                                                zoneNames);
  }
  else if (K_STRING::cmp(fileFmt, "fmt_stl") == 0) // fmt stl
  {
    if (fieldc.size() != 0)
      printf("Warning: convertArrays2File: structured arrays not converted.\n"); 
    
    isok = K_IO::GenIO::getInstance()->fstlwrite(fileName, dataFmt, varString,
                                                 ni, nj, nk,
                                                 fieldc, fieldu, connectu, elt,
                                                 zoneNames);
  }
  else if (K_STRING::cmp(fileFmt, "bin_3ds") == 0) // 3ds
  {
    if (fieldc.size() != 0)
      printf("Warning: convertArrays2File: structured arrays not converted.\n"); 
    
    isok = K_IO::GenIO::getInstance()->f3dswrite(fileName, dataFmt, varString,
                                                 ni, nj, nk,
                                                 fieldc, fieldu, connectu, elt,
                                                 zoneNames);
  }
  else if (K_STRING::cmp(fileFmt, "bin_ply") == 0) // ply
  {
    if (fieldc.size() != 0)
      printf("Warning: convertArrays2File: structured arrays not converted.\n"); 
    
    isok = K_IO::GenIO::getInstance()->plywrite(fileName, dataFmt, varString,
                                                ni, nj, nk,
                                                fieldc, fieldu, connectu, elt,
                                                zoneNames);
  }
  else if (K_STRING::cmp(fileFmt, "bin_wav") == 0) // bin wav
  { 
    isok = K_IO::GenIO::getInstance()->wavwrite(fileName, dataFmt, varString,
                                                ni, nj, nk,
                                                fieldc, fieldu, connectu, elt, 
                                                zoneNames);
  }
  else if (K_STRING::cmp(fileFmt, "fmt_xfig") == 0) // fmt xfig
  { 
    isok = K_IO::GenIO::getInstance()->xfigwrite(fileName, dataFmt, varString,
                                                 ni, nj, nk,
                                                 fieldc, fieldu, connectu, elt,
                                                 zoneNames);
  }
  else if (K_STRING::cmp(fileFmt, "fmt_svg") == 0) // fmt svg
  { 
    isok = K_IO::GenIO::getInstance()->svgwrite(fileName, dataFmt, varString,
                                                ni, nj, nk,
                                                fieldc, fieldu, connectu, elt,
                                                zoneNames);
  }
  else if (K_STRING::cmp(fileFmt, "fmt_cedre") == 0) // fmt cedre
  { 
    if (fieldc.size() != 0)
      printf("Warning: convertArrays2File: structured arrays not converted.\n"); 
    isok = K_IO::GenIO::getInstance()->cedrewrite(
      fileName, dataFmt, varString,
      ni, nj, nk,
      fieldc, fieldu, connectu, elt,
      zoneNames, BCFacesO);
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "convertArrays2File: unrecognised format.");
    return NULL;
  }
  if (isok == 1) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "convertArrays2File: file not written.");
    return NULL;
  }
  printf("done.\n"); fflush(stdout);
  
  // Deleting fields
  E_Int sizefieldc = fieldc.size();
  for (E_Int i = 0; i < sizefieldc; i++) delete fieldc[i];
  
  E_Int sizefieldu = fieldu.size();  
  for (E_Int i = 0; i < sizefieldu; i++)
  {
    delete fieldu[i]; delete connectu[i];
  }
  
  Py_INCREF(Py_None);
  return Py_None;
}

//=============================================================================
/* Check if file format is recognised
   Return 1 if file format is take into account
   Return 0 otherwise. */
//=============================================================================
E_Int K_CONVERTER::checkRecognisedFormat(char* fileFmt)
{
  if (K_STRING::cmp(fileFmt, "fmt_tp") == 0 ||
      K_STRING::cmp(fileFmt, "bin_tp") == 0 ||
      K_STRING::cmp(fileFmt, "fmt_v3d") == 0 ||
      K_STRING::cmp(fileFmt, "bin_v3d") == 0 ||
      K_STRING::cmp(fileFmt, "fmt_plot3d") == 0 ||
      K_STRING::cmp(fileFmt, "bin_plot3d") == 0 ||
      K_STRING::cmp(fileFmt, "fmt_pov") == 0 ||
      K_STRING::cmp(fileFmt, "bin_df3") == 0 ||
      K_STRING::cmp(fileFmt, "fmt_mesh") == 0 ||
      K_STRING::cmp(fileFmt, "fmt_gmsh") == 0 ||
      K_STRING::cmp(fileFmt, "bin_gmsh") == 0 ||
      K_STRING::cmp(fileFmt, "bin_stl") == 0 ||
      K_STRING::cmp(fileFmt, "fmt_stl") == 0 ||
      K_STRING::cmp(fileFmt, "bin_wav") == 0 ||
      K_STRING::cmp(fileFmt, "fmt_xfig") == 0 ||
      K_STRING::cmp(fileFmt, "fmt_svg") == 0 ||
      K_STRING::cmp(fileFmt, "fmt_obj") == 0 ||
      K_STRING::cmp(fileFmt, "fmt_gts") == 0 ||
      K_STRING::cmp(fileFmt, "fmt_iges") == 0 ||
      K_STRING::cmp(fileFmt, "bin_png") == 0 ||
      K_STRING::cmp(fileFmt, "bin_3ds") == 0 ||
      K_STRING::cmp(fileFmt, "bin_ply") == 0 ||
      K_STRING::cmp(fileFmt, "fmt_cedre") == 0 ||
      K_STRING::cmp(fileFmt, "fmt_tgf") == 0 ||
      K_STRING::cmp(fileFmt, "fmt_su2") == 0)
  {
    return 1;
  }
  return 0;
}

//=============================================================================
/* Returns the element type id */
//=============================================================================
E_Int K_CONVERTER::getElementTypeId(const char* eltType)
{
  if (K_STRING::cmp(eltType, "NODE") == 0) // NODE-> structure
    return 0;
  if (K_STRING::cmp(eltType, "BAR") == 0)
    return 1;
  if (K_STRING::cmp(eltType, "TRI") == 0)
    return 2;
  if (K_STRING::cmp(eltType, "QUAD") == 0)
    return 3;
  if (K_STRING::cmp(eltType, "TETRA") == 0)
    return 4;
  if (K_STRING::cmp(eltType, "PYRA") == 0)
    return 5;
  if (K_STRING::cmp(eltType, "PENTA") == 0)
    return 6;
  if (K_STRING::cmp(eltType, "HEXA") == 0)
    return 7;
  if (K_STRING::cmp(eltType, "NGON") == 0)
    return 8;
  return -1; //unknown
}

//=============================================================================
// Lit les paths specifies dans le fichier file (seult HDF).
// Retourne une liste d'objets pythons contenant les noeuds pointes par les
// chemins 
//=============================================================================
PyObject* K_CONVERTER::readPyTreeFromPaths(PyObject* self, PyObject* args)
{
  char* fileName; char* format; E_Int maxFloatSize; E_Int maxDepth; 
  PyObject* paths; PyObject* skipTypes;
  if (!PYPARSETUPLEI(args, "sOsllO", "sOsiiO", 
                     &fileName, &paths, &format, 
                     &maxFloatSize, &maxDepth, &skipTypes)) return NULL;
  
  if (skipTypes == Py_None) skipTypes = NULL;
  PyObject* ret = NULL;
  if (K_STRING::cmp(format, "bin_cgns") == 0)
    ret = K_IO::GenIO::getInstance()->hdfcgnsReadFromPaths(fileName, paths, maxFloatSize, maxDepth, skipTypes);
  else if (K_STRING::cmp(format, "bin_hdf") == 0)
    ret = K_IO::GenIO::getInstance()->hdfcgnsReadFromPaths(fileName, paths, maxFloatSize, maxDepth, skipTypes);
  else if (K_STRING::cmp(format, "bin_adf") == 0)
    ret = K_IO::GenIO::getInstance()->adfcgnsReadFromPaths(fileName, paths, maxFloatSize, maxDepth);
  else
  {
    PyErr_SetString(PyExc_TypeError, 
                    "readPyTreeFromPaths: unknown file format.");
    return NULL;
  }
  return ret;
} 

//=============================================================================
// Ecrit les paths specifies dans le fichier file (ADF/HDF).
//=============================================================================
PyObject* K_CONVERTER::writePyTreePaths(PyObject* self, PyObject* args)
{
  char* fileName; char* format; E_Int maxDepth; E_Int mode;
  PyObject* paths; PyObject* nodeList; PyObject* links;
  if (!PYPARSETUPLEI(args, "sOOsllO", "sOOsiiO", &fileName, &nodeList, &paths, &format, &maxDepth, &mode, &links))
    return NULL;
  
  if (links == Py_None) { links = NULL; }
  
  E_Int ret = 1;

  if (K_STRING::cmp(format, "bin_cgns") == 0)
    ret = K_IO::GenIO::getInstance()->hdfcgnsWritePaths(fileName, nodeList, paths, links, maxDepth, mode);
  else if (K_STRING::cmp(format, "bin_hdf") == 0)
    ret = K_IO::GenIO::getInstance()->hdfcgnsWritePaths(fileName, nodeList, paths, links, maxDepth, mode);
  else if (K_STRING::cmp(format, "bin_adf") == 0)
    ret = K_IO::GenIO::getInstance()->adfcgnsWritePaths(fileName, nodeList, paths, maxDepth, mode);
  if (ret == 1) return NULL; // exceptions deja levees

  Py_INCREF(Py_None);
  return Py_None;
} 

//=============================================================================
// Delete des paths du fichier (ADF/HDF).
//=============================================================================
PyObject* K_CONVERTER::deletePyTreePaths(PyObject* self, PyObject* args)
{
  char* fileName; char* format;
  PyObject* paths;
  if (!PYPARSETUPLEI(args, "sOs", "sOs", &fileName, &paths, &format))
    return NULL;
  
  E_Int ret = 1;
  if (K_STRING::cmp(format, "bin_cgns") == 0)
    ret = K_IO::GenIO::getInstance()->hdfcgnsDeletePaths(fileName, paths);
  else if (K_STRING::cmp(format, "bin_hdf") == 0)
    ret = K_IO::GenIO::getInstance()->hdfcgnsDeletePaths(fileName, paths);
  else if (K_STRING::cmp(format, "bin_adf") == 0)
    ret = K_IO::GenIO::getInstance()->adfcgnsDeletePaths(fileName, paths);
  if (ret == 1) return NULL; // exceptions deja levees

  Py_INCREF(Py_None);
  return Py_None;
}
