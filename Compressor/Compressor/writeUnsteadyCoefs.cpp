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
# include <stdio.h>
#include <map>
# include "compressor.h"
using namespace std;
using namespace K_FLD;

#define WRITE1(ptr_file,f,i1) do if (strcmp(f,"f") == 0) {fprintf(ptr_file,"%d\n",i1);} else {fwrite( &i1 , sizeof(int) , 1 , ptr_file);} while(0)
#define WRITE2(ptr_file,f,i1,i2) do {if (strcmp(f,"f") == 0) {fprintf(ptr_file,"%d %d\n",i1,i2);} else {fwrite( &i1 , sizeof(int) , 1 , ptr_file);fwrite( &i2 , sizeof(int) , 1 , ptr_file);}} while(0)
#define WRITE3(ptr_file,f,i1,i2,tabf)\
do\
{\
  if (strcmp(f,"f") == 0) {fprintf(ptr_file, "%d %d %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f\n",i1,i2,tabf[0],tabf[1],tabf[2],tabf[3],tabf[4],tabf[5],tabf[6]);}\
  else\
  {\
    fwrite( &i1 , sizeof(int) , 1 , ptr_file);fwrite( &i2 , sizeof(int) , 1 , ptr_file);\
    fwrite( &tabf , sizeof(E_Float) , 7 , ptr_file);\
  }\
}\
while(0)
#define WRITE3F(ptr_file,f,i1,i2,tabf,i3)\
do\
{\
  if (strcmp(f,"f") == 0) {fprintf(ptr_file, "%d %d %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %d\n",i1,i2,tabf[0],tabf[1],tabf[2],tabf[3],tabf[4],tabf[5],tabf[6],i3);}\
  else\
  {\
    fwrite( &i1 , sizeof(int) , 1 , ptr_file);fwrite( &i2 , sizeof(int) , 1 , ptr_file);\
    fwrite( &tabf , sizeof(E_Float) , 7 , ptr_file);\
    fwrite( &i3 , sizeof(int) , 1 , ptr_file);\
  }\
}\
while(0)
#define WRITE4(ptr_file,f,i1,i2,i3,tabf)\
do\
{\
  if (strcmp(f,"f") == 0) {fprintf(ptr_file, "%d %d %d %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f\n",i1,i2,i3,tabf[0],tabf[1],tabf[2],tabf[3],tabf[4],tabf[5],tabf[6]);}\
  else\
  {\
    fwrite( &i1 , sizeof(int) , 1 , ptr_file);\
    fwrite( &i2 , sizeof(int) , 1 , ptr_file);\
    fwrite( &i3 , sizeof(int) , 1 , ptr_file);\
    fwrite( &tabf , sizeof(E_Float) , 7 , ptr_file);\
  }\
}\
while(0)
#define WRITE4F(ptr_file,f,i1,i2,i3,tabf,i4)\
do\
{\
  if (strcmp(f,"f") == 0) {fprintf(ptr_file, "%d %d %d %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %d\n",i1,i2,i3,tabf[0],tabf[1],tabf[2],tabf[3],tabf[4],tabf[5],tabf[6],i4);}\
  else\
  {\
    fwrite( &i1 , sizeof(int) , 1 , ptr_file);\
    fwrite( &i2 , sizeof(int) , 1 , ptr_file);\
    fwrite( &i3 , sizeof(int) , 1 , ptr_file);\
    fwrite( &tabf , sizeof(E_Float) , 7 , ptr_file);\
    fwrite( &i4 , sizeof(int) , 1 , ptr_file);\
  }\
}\
while(0)
#define WRITE5(ptr_file,f,i1,i2,i3,i4,tabf)\
do\
{\
  if (strcmp(f,"f") == 0) {fprintf(ptr_file, "%d %d %d %d %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f\n",i1,i2,i3,i4,tabf[0],tabf[1],tabf[2],tabf[3],tabf[4],tabf[5],tabf[6]);}\
  else\
  {\
    fwrite( &i1 , sizeof(int) , 1 , ptr_file);\
    fwrite( &i2 , sizeof(int) , 1 , ptr_file);\
    fwrite( &i3 , sizeof(int) , 1 , ptr_file);\
    fwrite( &i4 , sizeof(int) , 1 , ptr_file);\
    fwrite( &tabf , sizeof(E_Float) , 7 , ptr_file);\
  }\
}\
while(0)
#define WRITE5F(ptr_file,f,i1,i2,i3,i4,tabf,i5)\
do\
{\
  if (strcmp(f,"f") == 0) {fprintf(ptr_file, "%d %d %d %d %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %d\n",i1,i2,i3,i4,tabf[0],tabf[1],tabf[2],tabf[3],tabf[4],tabf[5],tabf[6],i5);}\
  else\
  {\
    fwrite( &i1 , sizeof(int) , 1 , ptr_file);\
    fwrite( &i2 , sizeof(int) , 1 , ptr_file);\
    fwrite( &i3 , sizeof(int) , 1 , ptr_file);\
    fwrite( &i4 , sizeof(int) , 1 , ptr_file);\
    fwrite( &tabf , sizeof(E_Float) , 7 , ptr_file);\
    fwrite( &i5 , sizeof(int) , 1 , ptr_file);\
  }\
}\
while(0)

//=============================================================================
/* Ecriture des coefficients d'interpolation dans un fichier relisible 
   par elsA dans un format optimise pour l instationnaire */
//=============================================================================
PyObject* K_COMPRESSOR::writeUnsteadyCoefs(PyObject* self, PyObject* args)
{
  IMPORTNUMPY;
  E_Int Iteration;          // tableau contenant le numero des iterations
  PyObject *ListInterpData; // donnees d interpolations pour un donneur
  char* FileName;           // nom du fichier
  char* Localisation;       // cellules ou centres des faces
  char* Format;             // format du fichier

  if (!PYPARSETUPLEI(args, "lOsss", "iOsss",
                        &Iteration, &ListInterpData, &FileName, &Localisation, &Format)) 
    return NULL;
  /*-------------------*/
  /* Variables locales */
  /*-------------------*/
  E_Int intKey1, intKey2, nbPts;

  /*----------------------------------------*/
  /* Extraction de l iteration */
  /*----------------------------------------*/
  E_Int iteration = Iteration;

  /*----------------------------------------*/
  /* Extraction des donnees d interpolation */
  /*----------------------------------------*/
  list<map<E_Int,map<E_Int,vector<E_Float> > > > listInterpData;
  if (PyList_Check (ListInterpData) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "deltaInterpolations: listInterpData must be a list.");
    return NULL;    
  }
  E_Int nbIterations = PyList_Size(ListInterpData);
  for (E_Int ite  = 0 ; ite < nbIterations; ite++)
  {
    PyObject* InterpData = PyList_GetItem(ListInterpData, ite);
    if (PyDict_Check (InterpData) == 0)
    {
      PyErr_SetString(PyExc_TypeError, 
                      "deltaInterpolations: listInterpData must be a list of dictionary.");
      return NULL;
    }
    PyObject * key_list = PyDict_Keys(InterpData);
    E_Int nbRcvZones = PyList_Size(key_list);
    //Get keys and corresponding values from the dictionary
    map<E_Int,map<E_Int, vector<E_Float> > > mapInterpDataByBlock;
    for (E_Int i  = 0 ; i < nbRcvZones; i++)
    {
      PyObject * pyKey = PyList_GetItem(key_list, i);
      if (PyInt_Check(pyKey) == 0)
      {
        printf("Warning: deltaInterpolations: invalid int for variable %d. Skipped...\n", i);
      }
      else
        intKey1 = PyInt_AsLong(pyKey);
      // map <E_Int, <listE_Float > >
      PyObject* InterpDataByRcvBlk = PyDict_GetItem(InterpData, pyKey);
      if (PyDict_Check (InterpDataByRcvBlk) == 0)
      {
        PyErr_SetString(PyExc_TypeError, 
                        "deltaInterpolations: listInterpData must be a list of dictionary.");
        return NULL;
      }
      PyObject * key_list2 = PyDict_Keys(InterpDataByRcvBlk);
      E_Int nbIndices = PyList_Size(key_list2);
      //DBG TO COMMENT
      map<E_Int, vector<E_Float> > mapInterpDataByIndex;
      for (E_Int iterIndex  = 0 ; iterIndex < nbIndices; iterIndex++)
      {
        PyObject * keyIndexRcv = PyList_GetItem(key_list2, iterIndex);
        if (PyInt_Check(keyIndexRcv) == 0)
        {
          printf("Warning: deltaInterpolations: invalid int for variable %d. Skipped ...\n", iterIndex);
        }
        else
          intKey2 = PyInt_AsLong(keyIndexRcv);
        // vector<E_Float>
        vector<E_Float> listOfData;
        PyObject* pyListValues = PyDict_GetItem(InterpDataByRcvBlk, keyIndexRcv);
        // check if pyListValue is a list
        if (PyList_Check (pyListValues) == 0)
        {
          PyErr_SetString(PyExc_TypeError, 
                          "deltaInterpolations: listInterpData must be a list of dictionary.");
          return NULL;
        }
        // list<E_Float>
        E_Int tmpListValueSize = PyList_Size(pyListValues);
        for (E_Int v  = 0 ; v < tmpListValueSize; v++)
        {
          PyObject* pyValue = PyList_GetItem(pyListValues, v);
          // // DBG
          // PyTypeObject* type = pyValue->ob_type;
          // const char* p = type->tp_name;
          // printf("DBG pyValue %s\n",p);fflush(stdout);
          // DBG
          if (PyInt_Check(pyValue) != 0)
          {
            E_Int intValue = PyInt_AsLong(pyValue);
            listOfData.push_back((E_Float)(intValue));
            // DBG
            //printf("DBG value = %d\n",intValue);fflush(stdout);
            // DBG
          }
          // DBG
          // else if (PyArray_IsScalar(pyValue, Integer))
          // {
          //   E_Int intValue = PyInt_AsLong(pyValue);
          //   listOfData.push_back((E_Float)(intValue));
          // }
          // DBG
          else if (PyFloat_Check(pyValue) != 0)
          {
            E_Float floatValue = PyFloat_AsDouble(pyValue);
            listOfData.push_back(floatValue);
          }
          else
          {
            printf("deltaInterpolations: stored data [storage flag, donor index, periodicity, coefficients] must be integer of float values.\n");
          }
        }
        mapInterpDataByIndex[intKey2] = listOfData;
      }
      mapInterpDataByBlock[intKey1]=mapInterpDataByIndex;
    }
    listInterpData.push_back(mapInterpDataByBlock);
  }
 
  printf("DBG Ecriture des fichiers d interpolation\n");fflush(stdout);
  /*--------------------------------------------------------*/
  /* Ecriture des fichiers d interpolation                  */
  /*--------------------------------------------------------*/
  // Ouverture du fichier
  FILE* ptr_file = NULL;
  char* file = new char[K_ARRAY::VARSTRINGLENGTH];
  strcpy(file,FileName);
  char* loc = new char[K_ARRAY::VARSTRINGLENGTH];
  strcpy(loc,Localisation);
  char* format = new char[K_ARRAY::VARSTRINGLENGTH];
  strcpy(format,Format);
  ptr_file = fopen(file, "w");
  // Ecriture de la premiere iteration de lecture
  printf("DBG iteration %d\n",iteration);fflush(stdout);
  printf("DBG loc %s\n",loc);fflush(stdout);
  printf("DBG format = %s, test = %d\n",format,strcmp(format,"f"));fflush(stdout);
  WRITE1(ptr_file,format,iteration);
  
    // Parcours des iterations
    for (list<map<E_Int, map <E_Int, vector<E_Float> > > >::iterator itr = listInterpData.begin(); itr!= listInterpData.end(); itr++)
    {
      map<E_Int, map <E_Int, vector<E_Float> > > interpData = *itr;
      // Ecriture du nombre de blocs receveur pour l iteration courante
      E_Int nbRcvZones = interpData.size();
      WRITE1(ptr_file,format,nbRcvZones);
      // Parcours des zones receveuses
      for (map<E_Int, map <E_Int, vector<E_Float> > >::iterator iterz = interpData.begin();
           iterz != interpData.end();
           iterz++)
      {
        E_Int rcvId = (*iterz).first;
        WRITE1(ptr_file,format,rcvId);
        map <E_Int, vector<E_Float> > interpDataByRcvBlk = (*iterz).second;
        // Ecriture du nombre de points interpoles
        nbPts = interpDataByRcvBlk.size();
        WRITE1(ptr_file,format,nbPts);
        E_Int indexDonor;
        E_Int periodicType; 
        E_Int faceDir;
        FldArrayF interpCoefs(7); // DBG FAIRE PLUS GENERAL
        // Parcours des donnees d interpolation par bloc
        for (map <E_Int, vector<E_Float> >::iterator itri  = interpDataByRcvBlk.begin();
             itri != interpDataByRcvBlk.end(); 
             itri++)
        {
          // indice receveur
          E_Int indexRcv = (*itri).first;
          vector<E_Float>  interpDataByIndex = (*itri).second;
          // Flag determinant les donnees stockees
          E_Int flag = interpDataByIndex[0];
          // direction de la face par rapport a la cellule frontiere adjacente
          if (!strcmp(loc,"face")) faceDir = (E_Int)(interpDataByIndex[interpDataByIndex.size()-1]);
        
          switch (flag)
          {
            case 0: // nouveau point ou point avec un donneur different (stockage complet)
              indexDonor = (E_Int)(interpDataByIndex[1]);
              periodicType = (E_Int)(interpDataByIndex[2]);
              for (E_Int cf=0; cf < 7; cf++)
                interpCoefs[cf] = interpDataByIndex[cf+3];
              if (!strcmp(loc,"cell")) WRITE5(ptr_file,format,flag, indexRcv, indexDonor, periodicType, interpCoefs);
              else WRITE5F(ptr_file,format,flag, indexRcv, indexDonor, periodicType, interpCoefs, faceDir);
                break;
            case -1: // point a detruire
              WRITE2(ptr_file,format,flag, indexRcv);
              break;
            case 1: // meme point receveur, meme donneur, meme periodicite
              for (E_Int cf=0; cf < 7; cf++)
                interpCoefs[cf] = interpDataByIndex[cf+1];
              if (!strcmp(loc,"cell")) WRITE3(ptr_file,format,flag, indexRcv, interpCoefs);
              else WRITE3F(ptr_file,format,flag, indexRcv, interpCoefs, faceDir);
                break;
            case 2: // meme point receveur, meme donneur, la periodicite a change
              periodicType = (E_Int)(interpDataByIndex[1]);
              for (E_Int cf=0; cf < 7; cf++)
                interpCoefs[cf] = interpDataByIndex[cf+2];
              if (!strcmp(loc,"cell"))  WRITE4(ptr_file,format,flag, indexRcv, periodicType, interpCoefs);
              else WRITE4F(ptr_file,format,flag, indexRcv, periodicType, interpCoefs, faceDir);
                break;
            case 3: // donneur different, la periodicite est la meme
              indexDonor = (E_Int)(interpDataByIndex[1]);
              for (E_Int cf=0; cf < 7; cf++)
                interpCoefs[cf] = interpDataByIndex[cf+2];
              if (!strcmp(loc,"cell"))  WRITE4(ptr_file,format,flag, indexRcv, indexDonor, interpCoefs);
              else  WRITE4F(ptr_file,format,flag, indexRcv, indexDonor, interpCoefs, faceDir);
                break;
            case 4: // memes donnees
              WRITE2(ptr_file,format,flag, indexRcv);
                break;
            case 5: // meme point receveur, meme donneur, meme periodicite. Only for face
              WRITE3F(ptr_file,format,flag, indexRcv, interpCoefs, faceDir);
                break;
            case 6: // meme point receveur, meme donneur. Only for face
              WRITE4F(ptr_file,format,flag, indexRcv, periodicType, interpCoefs, faceDir);
                break;
            default:
              PyErr_SetString(PyExc_TypeError,"deltaInterpolations: bad flag for interpolation storage.");
              break;
          } // Fin switch flag
        } // Fin de la boucle des donnees d interpolation par bloc
      } // Fin de la boucle sur les zones receveuses
    } // End loop on iterations

  // Fermeture du fichier
  fclose(ptr_file);

  Py_INCREF(Py_None);
  return Py_None;
}
