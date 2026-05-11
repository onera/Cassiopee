/*
    Copyright 2013-2025 Onera.

    This file is part of the FFX (Far-Field Exergy) post-processing module.
*/
# include <string.h>
# include <math.h>

using namespace K_FLD;
using namespace std;

//=============================================================================
/* Add layers around a control volume defined by a tag variable.
   The tag for the volume is located on cell centers. */
//=============================================================================
PyObject* K_FFX::addLayers(PyObject* self,PyObject* args)
{
  PyObject* array; char* tagname;
  E_Int nlayers;
  if (!PYPARSETUPLE_(args, O_ S_ I_,
                    &array, &tagname, &nlayers)) return NULL;

  // Check array
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cn;
  E_Int ni, nj, nk;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk,
                                     cn, eltType);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "addLayers: invalid array.");
    return NULL;
  }

  E_Int postag = K_ARRAY::isNamePresent(tagname, varString);
  if (postag == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "addLayers: tag variable not found in array.");
    RELEASESHAREDB(res,array,f,cn); return NULL;
  }
  postag++;

  // Add extra layers aroung the volume defined by the tag variable
  if (res == 1) // struct
  {
    //E_Int nic = K_FUNC::E_max(1,ni-1);
    //E_Int njc = K_FUNC::E_max(1,nj-1);
    //E_Int nkc = K_FUNC::E_max(1,nk-1);
    //E_Int ncells = nic*njc*nkc;
    //E_Int nicnjc = nic*njc;
    //E_Int icc;
    //for (E_Int k = 0; k < nkc; k++)
    //  for (E_Int j = 0; j < njc; j++)
    //    for (E_Int i = 1; i < nic; i++)
    //    {
    //      icc = i+j*ni+k*ninjc;
    //    }
    PyErr_SetString(PyExc_TypeError,
                    "addLayers: only implemented for NGON unstructured zones.");
    RELEASESHAREDB(res,array,f,cn); return NULL;
  }
  else if (res == 2) // unstruct
  {
    //if (strcmp(eltType, "NGON") == 0 || strcmp(eltType, "NGON*") == 0)
    //if (strcmp(eltType, "NGON") == 0)
    if (strcmp(eltType, "NGON*") == 0)
    {
      FldArrayI cFE;
      K_CONNECT::connectNG2FE(*cn, cFE); // get face/elmt connectivity
      E_Int* cFE1 = cFE.begin(1);
      E_Int* cFE2 = cFE.begin(2);

      E_Int nfaces = cn->getNFaces();
      E_Int nelts = cn->getNElts();
      E_Float* ftp = f->begin(postag);

      E_Int e1, e2;
      for (E_Int l = 0; l < nlayers; l++)
      {
        for (E_Int i = 0; i < nfaces; i++)
        {
          e1 = cFE1[i]-1; e2 = cFE2[i]-1;
          if (e1 != -1 && e2 != -1) // if face has 2 neighboring elmts
          {
            if (ftp[e1] == K_CONST::ONE && ftp[e2] == K_CONST::E_ZERO_FLOAT)
              ftp[e2] = K_CONST::TWO;
            else if (ftp[e2] == K_CONST::ONE && ftp[e1] == K_CONST::E_ZERO_FLOAT)
              ftp[e1] = K_CONST::TWO;
          }
        }
        for (E_Int i = 0; i < nelts; i++)
        {
          if (ftp[i] == K_CONST::TWO) ftp[i] = K_CONST::ONE;
        }
      }
    }
    else
    {
      if (strcmp(eltType, "BAR") != 0 &&
          strcmp(eltType, "TRI") != 0 &&
          strcmp(eltType, "QUAD") != 0 &&
          strcmp(eltType, "TETRA") != 0 &&
          strcmp(eltType, "HEXA") != 0 &&
          strcmp(eltType, "PYRA") != 0 &&
          strcmp(eltType, "PENTA") != 0 )
      {
        PyErr_SetString(PyExc_TypeError,
                        "addLayers: not a valid element type.");
        RELEASESHAREDU(array,f,cn); return NULL;
      }
      PyErr_SetString(PyExc_TypeError,
                      "addLayers: only implemented for NGON unstructured zones.");
      RELEASESHAREDB(res,array,f,cn); return NULL;
    }
  }
  RELEASESHAREDB(res,array,f,cn)
  Py_INCREF(Py_None);
  return Py_None;
}
