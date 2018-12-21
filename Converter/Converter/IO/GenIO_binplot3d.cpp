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
// nasa plot3d binary file support

# include "Array/Array.h"
# include "GenIO.h"
# include <stdio.h>
# include <string.h>

using namespace K_ARRAY;
using namespace K_FLD;
using namespace std;

//=============================================================================
/* 
   Find xyz, q, and f file name if they exists in file.
   Its the caller responsability to free the returned char*.
   IN: file: file to be opened
   OUT: xyzFile, qFile, fFile: real file name that can be opened
   Retourne le nombre de fichiers que l'on peut ouvrir (0,1,2 ou 3).
   
   Regle:
   Si le fichier file existe, on determine si c'est un xyz, q ou f file
   suivant l'extension:
   - .gbin: xyz
   - .qbin: q
   - .fbin: f
   Sinon, on compose un nom file + .gbin, file + .qbin, file + .fbin, et
   on regarde si ces fichiers existent.
*/
//=============================================================================
E_Int K_IO::GenIO::findOpenablePlot3dFiles(
  char* file,
  char*& xyzFile, char*& qFile, char*& fFile)
{
  FILE* ptrFile;
  xyzFile = NULL; qFile = NULL; fFile = NULL;
  char temp[120];
  E_Boolean used;
  
  E_Int openable = 0;
  
  // Analyse la chaine file, decoupe suivant les virgules
  char split[3][120];
  E_Int c = 0;
  E_Int n = 0;
  E_Int l = 0;
  while (file[c] != '\0')
  {
    split[n][l] = file[c];
    if (file[c] == ',') // multiple files in string
    {
      split[n][l] = '\0';
      n++;
      if (n >= 4)
      {
        printf("Warning: too many file names provided. Only 3 are taken into account.\n");
        n--;
        break;
      }
      l = 0;
    }
    else
      l++;
    c++;
  }
  split[n][l] = '\0';
  
  // Examine chaque nom de fichier, existe-t-il?
  for (E_Int i = 0; i <= n; i++)
  {
    used = false;
    l = strlen(split[i]);
    
    ptrFile = fopen(split[i], "rb");

    // Le fichier peut-etre ouvert
    if (ptrFile != NULL)
    {
      used = true;
      if (l >= 4)
      {
        temp[0] = split[i][l-4]; temp[1] = split[i][l-3];  
        temp[2] = split[i][l-2]; temp[3] = split[i][l-1];
        temp[4] = '\0';
        if (strcmp(temp, "gbin") == 0)
        {
          if (xyzFile == NULL)
          {
            xyzFile = new char[120];
            openable++;
          }
          else
            printf("Warning: multiple file names given for xyz definition. Discarding previous one.\n");
          strcpy(xyzFile, split[i]);
        }
        else if (strcmp(temp, "qbin") == 0)
        {
          if (qFile == NULL)
          {
            qFile = new char[120];
            openable++;
          }
          else
            printf("Warning: multiple file names given for q definition. Discarding previous one.\n");
          strcpy(qFile, split[i]);
        }
        else if (strcmp(temp, "fbin") == 0)
        {
          if (fFile == NULL)
          {
            fFile = new char[120];
            openable++;
          }
          else
            printf("Warning: multiple file names given for f definition. Discarding previous one.\n");
          strcpy(fFile, split[i]);
        }
        else
        {
          printf("Warning: findOpenablePlot3dFiles: assuming that file is an xyz plot 3d file.\n");
          xyzFile = new char[120];
          openable++;
          strcpy(xyzFile, split[i]);
        }
      }
      else
      {
        printf("Warning: findOpenablePlot3dFiles: assuming that file is an xyz  plot 3d file.\n");
        xyzFile = new char[120];
        openable++;
        strcpy(xyzFile, split[i]);
      }
      fclose(ptrFile);
    }
    else
    {
      // Le nom n'existe pas. On verifie s'il s'agit d'un nom generique
      strcpy(temp, split[i]);
      strcat(temp, ".gbin");
      ptrFile = fopen(temp, "rb");
      
      if (ptrFile != NULL)
      {
        used = true;
        xyzFile = new char[120];
        openable++;
        strcpy(xyzFile, temp);
        fclose(ptrFile);
      }
      
      strcpy(temp, split[i]);
      strcat(temp, ".qbin");
      ptrFile = fopen(temp, "rb");
      
      if (ptrFile != NULL)
      {
        used = true;
        qFile = new char[120];
        openable++;
        strcpy(qFile, temp);
        fclose(ptrFile);
      }
      
      strcpy(temp, split[i]);
      strcat(temp, ".fbin");
      ptrFile = fopen(temp, "rb");
      
      if (ptrFile != NULL)
      {
        used = true;
        fFile = new char[120];
        openable++;
        strcpy(fFile, temp);
        fclose(ptrFile);
      }
    }
    if (used == false)
      printf("Warning: findOpenablePlot3dFiles: the name %s in file string can not be used.\n", split[i]);
  } // pour chaque nom splitte
  
  return openable;
}

//=============================================================================
// Trouve les noms reels de fichiers a ecrire.
// IN: file: nom du fichier
// IN: varString: chaine des variables
// OUT: xyzFile, qFile, fFile: nom des fichiers reels a ecrire.
//=============================================================================
E_Int K_IO::GenIO::findWritablePlot3dFiles(
  char* file,
  char* varString,
  char*& xyzFile, char*& qFile, char*& fFile)
{
  char temp[120];
  E_Int writable = 0;

  // Analyse la chaine file, decoupe suivant les virgules
  char split[3][120];
  E_Int c = 0;
  E_Int n = 0;
  E_Int l = 0;
  while (file[c] != '\0')
  {
    split[n][l] = file[c];
    if (file[c] == ',') // multiple files in string
    {
      split[n][l] = '\0';
      n++;
      if (n >= 4)
      {
        printf("Warning: you have provided too much file names. Only 3 are taken into account.\n");
        n--;
        break;
      }
      l = 0;
    }
    else
      l++;
    c++;
  }
  split[n][l] = '\0';
  
  E_Int resx = isCoordinateXPresent(varString);
  E_Int resy = isCoordinateYPresent(varString);
  E_Int resz = isCoordinateZPresent(varString);
  //E_Int rescelln = isNamePresent("cellN", varString);
  //if (rescelln == -1) rescelln = isNamePresent("cellnf", varString);
  E_Int resro = isDensityPresent(varString);
  E_Int resrou = isMomentumXPresent(varString);
  E_Int resrov = isMomentumYPresent(varString);
  E_Int resrow = isMomentumZPresent(varString);
  E_Int resroE = isEnergyStagnationDensityPresent(varString);

  // Examine le nom du fichier
  for (E_Int i = 0; i <= n; i++)
  {
    // xyz files
    if (resx != -1 && resy != -1 && resz != -1)
    {
      l = strlen(split[i]);
      if (l >= 4)
      {
        temp[0] = split[i][l-4]; temp[1] = split[i][l-3];  
        temp[2] = split[i][l-2]; temp[3] = split[i][l-1];
        temp[4] = '\0';
        if (strcmp(temp, "gbin") == 0)
        {
          if (xyzFile == NULL)
          {
            xyzFile = new char[120];
            writable++;
          }
          else
            printf("Warning: you gave multiple file names for xyz definition. Discarding previous one.\n");
          strcpy(xyzFile, split[i]);
        }
      }
    }

    // q file
    if (resro != -1 && resrou != -1 && resrov != -1 &&
        resrow != -1 && resroE != -1 )
    {
      l = strlen(split[i]);
      if (l >= 4)
      {
        temp[0] = split[i][l-4]; temp[1] = split[i][l-3];  
        temp[2] = split[i][l-2]; temp[3] = split[i][l-1];
        temp[4] = '\0';
        if (strcmp(temp, "qbin") == 0)
        {
          if (qFile == NULL)
          {
            qFile = new char[120];
            writable++;
          }
          else
            printf("Warning: you gave multiple file names for q definition. Discarding previous one.\n");
          strcpy(qFile, split[i]);
        }
      }
    }
    
    // f file
    l = strlen(split[i]);
    if (l >= 4)
    {
      temp[0] = split[i][l-4]; temp[1] = split[i][l-3];  
      temp[2] = split[i][l-2]; temp[3] = split[i][l-1];
      temp[4] = '\0';
      if (strcmp(temp, "fbin") == 0)
      {
        if (fFile == NULL)
        {
          fFile = new char[120];
          writable++;
        }
        else
          printf("Warning: you gave multiple file names for f definition. Discarding previous one.\n");
        strcpy(fFile, split[i]);
      }
    }
  } // pour chaque nom splitte
  
  // Examine le nom du fichier
  for (E_Int i = 0; i <= n; i++)
  {
    if (xyzFile == NULL && resx != -1 && resy != -1 && resz != -1)
    {
      // nom generique
      xyzFile = new char[120];
      strcpy(xyzFile, split[i]);
      strcat(xyzFile, ".gbin");
      writable++;
    }
    
    if (qFile == NULL && resro != -1 && resrou != -1 && resrov != -1 &&
        resrow != -1 && resroE != -1 )
    {
      qFile = new char[120];
      writable++;
      strcpy(qFile, split[i]);
      strcat(qFile, ".qbin");
    }

    if (fFile == NULL)
    {
      fFile = new char[120];
      writable++;
      strcpy(fFile, split[i]);
      strcat(fFile, ".fbin");
    }
  }
  
  return writable;
}

//=============================================================================
/* 
   Read xyz file.
 */
//=============================================================================
E_Int K_IO::GenIO::xyzread(
  char* file, char*& varString,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& field)
{
  FILE* ptrFile;
  E_Int s;
  int ib;

  /* Constants */
  E_Int si = sizeof(int);
  E_Int si2 = sizeof(E_LONG);
  
  // Check endianess
  E_Int ret = v3dCheckEndian(file); // similar to v3d
  if (ret == -1)
  {
    printf("Warning: xyzread: can not open file %s.\n", file);
    return 1;
  }
  if (ret == -2)
  {
    printf("Warning: xyzread: format in file %s is not recognised.\n", file);
    return 1;
  }
  else if (ret == 0)
  {
    _convertEndian = false;
    _intLength = 4;
  }
  else if (ret == 1)
  {
    _convertEndian = false;
    _intLength = 8;
  }
  else if (ret == 2)
  {
    _convertEndian = true;
    _intLength = 4;
  }
  else
  {
    _convertEndian = true;
    _intLength = 8;
  }

  /* Opening */
  ptrFile = fopen(file, "rb");

  // Read number of zones in file
  fread(&ib, si, 1, ptrFile); // sep : taille int
  E_Int nzones;
  readInt(nzones, _intLength, 
          ptrFile, _convertEndian, si, si2);
  readInt(s, si, // sep : taille int
          ptrFile, _convertEndian, si, si2);
  if (s != _intLength)
  {
    printf("Warning: binary file corrupted.\n");
    fclose(ptrFile);
    return 1;
  }

  readInt(s, si, 
          ptrFile, _convertEndian, si, si2);
  if (s != _intLength*3*nzones) // sep : nombre de zones * 3 * taille int
  {
    printf("Warning : binary file corrupted.\n");
    fclose(ptrFile);
    return 1;
  }
  
  // Read dimensions of zones
  for (E_Int i = 0; i < nzones; i++)
  {
    E_Int ip, jp, kp;
    readInt(ip, _intLength, 
            ptrFile, _convertEndian, si, si2);    
    readInt(jp, _intLength, 
            ptrFile, _convertEndian, si, si2);    
    readInt(kp, _intLength, 
            ptrFile, _convertEndian, si, si2);

    ni.push_back(ip);
    nj.push_back(jp);
    nk.push_back(kp);
  }
  
  E_Int iblank = false;

  // Sep
  readInt(s, si, // sep : nombre de zone * taille int * 3
          ptrFile, _convertEndian, si, si2);    
  if (s != _intLength*3*nzones)
  {
    printf("Warning: binary file corrupted.\n");
    fclose(ptrFile);
    return 1;
  }

  // Read x,y,z coordinates
  for (E_Int i = 0; i < nzones; i++)
  {
    readInt(s, si, // sep taille dom * 3 * taille float
          ptrFile, _convertEndian, si, si2);
    E_Int np = ni[i]*nj[i]*nk[i];

    if (s == np*(3*8+_intLength))
    {
      _realLength = 8;
      iblank = true;
    }
    else if (s == np*(3*4+_intLength))
    {
      _realLength = 4;
      iblank = true;
    }
    else if (s == np*3*8)
    {
      _realLength = 8;
      iblank = false;
    }
    else if (s == np*3*4)
    {
      _realLength = 4;
      iblank = false;
    }
    else
    {
      printf("Warning: the file is corrupted.\n");
      fclose(ptrFile);
      return 1;
    }
   
    E_Int ip = ni[i];
    E_Int jp = nj[i];
    E_Int kp = nk[i];
    E_Int ipjpkp = ip*jp*kp;
    FldArrayF* f = NULL;
    if (iblank == true)
    {
      f = new FldArrayF(ipjpkp, 4);
      if (_realLength == 8)
        fread(f->begin(), sizeof(E_Float), 3*ipjpkp, ptrFile);
      else
      {
        float* buf = new float[3*ipjpkp];
        fread(buf, sizeof(float), 3*ipjpkp, ptrFile);
        E_Float* pf = f->begin();
        for (E_Int ind = 0; ind < 3*ipjpkp; ind++) pf[ind] = buf[ind];
        delete [] buf;
      }
      if (_intLength == 8)
      {
        E_LONG* blank = new E_LONG[ipjpkp];
        fread(blank, sizeof(E_LONG), ipjpkp, ptrFile);
        E_Float* pf = f->begin();
        for (E_Int ind = 0; ind < ipjpkp; ind++)
          pf[3*ipjpkp+ind] = blank[ind];
        delete [] blank;
      }
      else
      {
        int* blank = new int[ipjpkp];
        fread(blank, sizeof(int), ipjpkp, ptrFile);
        E_Float* pf = f->begin();
        for (E_Int ind = 0; ind < ipjpkp; ind++)
          pf[3*ipjpkp+ind] = blank[ind];
        delete [] blank;
      }
    }
    else
    {
      f = new FldArrayF(ipjpkp, 3);
      fread(f->begin(), sizeof(E_Float), 3*ipjpkp, ptrFile);
    }
    field.push_back(f);
    
    readInt(s, si, // sep taille dom *3 ou *4 taille float
          ptrFile, _convertEndian, si, si2);
  }

  fclose(ptrFile);

  // Build varString
  if (iblank == true)
    strcpy(varString, "x,y,z,cellN");
  else
    strcpy(varString, "x,y,z");
  return 0;
}

//=============================================================================
/* Read q file. */
//=============================================================================
E_Int K_IO::GenIO::p3dqread(
  char* file, char*& varString,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& field)
{
  FILE* ptrFile;
  E_Int s;
  int ib;

  E_Int xyzdone = false;
  if (ni.size() != 0) xyzdone = true;

  /* Constants */
  E_Int si = sizeof(int);
  E_Int si2 = sizeof(E_LONG);
  
  // Check endianess
  E_Int ret = v3dCheckEndian(file); // similar to v3d
  if (ret == -1)
  {
    printf("Warning: qread: can not open file %s.\n", file);
    return 1;
  }
  if (ret == -2)
  {
    printf("Warning: qread: format in file %s is not recognised.\n", file);
    return 1;
  }
  else if (ret == 0)
  {
    _convertEndian = false;
    _intLength = 4;
  }
  else if (ret == 1)
  {
    _convertEndian = false;
    _intLength = 8;
  }
  else if (ret == 2)
  {
    _convertEndian = true;
    _intLength = 4;
  }
  else
  {
    _convertEndian = true;
    _intLength = 8;
  }

  /* Opening */
  ptrFile = fopen(file, "rb");

  // Read number of zones in file
  fread(&ib, si, 1, ptrFile); // sep
  E_Int nzones;
  readInt(nzones, _intLength, 
          ptrFile, _convertEndian, si, si2);
  readInt(s, si, // sep : taille int
          ptrFile, _convertEndian, si, si2);
  if (s != _intLength)
  {
    printf("Warning: binary file corrupted.\n");
    fclose(ptrFile);
    return 1;
  }

  // Coherence
  if (xyzdone == true)
  {
    E_Int nisize = ni.size();
    if (nisize != nzones)
      printf("Warning: qread: the number of zones in q and xyz files is different.\n");
  }

  readInt(s, si, 
          ptrFile, _convertEndian, si, si2);
  if (s != _intLength*3*nzones) // sep : nombre de zones * 3 * taille int
  {
    printf("Warning: binary file corrupted.\n");
    fclose(ptrFile);
    return 1;
  }
  
  // Read dimensions of zones
  for (E_Int i = 0; i < nzones; i++)
  {
    E_Int ip, jp, kp;
    readInt(ip, _intLength, 
            ptrFile, _convertEndian, si, si2);    
    readInt(jp, _intLength, 
            ptrFile, _convertEndian, si, si2);    
    readInt(kp, _intLength, 
            ptrFile, _convertEndian, si, si2);

    // Coherence
    if (xyzdone == true)
    {
      if (ip != ni[i] || jp != nj[i] || kp != nk[i])
        printf("Warning: zones in q file and xyz are different.\n");
    }
    else
    {
      ni.push_back(ip);
      nj.push_back(jp);
      nk.push_back(kp);
    }
  }
  
  // Sep
  readInt(s, si, 
          ptrFile, _convertEndian, si, si2);
  if (s != _intLength*3*nzones) // sep : nombre de zones * 3 * taille int
  {
    printf("Warning: binary file corrupted.\n");
    fclose(ptrFile);
    return 1;
  }

  // Read ro, rou, rov, row, roE
  for (E_Int i = 0; i < nzones; i++)
  {
    readInt(s, si, 
            ptrFile, _convertEndian, si, si2);
    if (s != 4*sizeof(E_Float)) // sep 4 * taille float
    {
      printf("Warning: binary file corrupted.\n");
      fclose(ptrFile);
      return 1;
    }

    for (E_Int p = 0; p < 8; p++)
      fread(&ib, si, 1, ptrFile); // Mach, etc..
    
    readInt(s, si, 
            ptrFile, _convertEndian, si, si2);
    if (s != 4*sizeof(E_Float)) // sep 4 * taille float
    {
      printf("Warning: binary file corrupted.\n");
      fclose(ptrFile);
      return 1;
    }

    E_Int ip = ni[i];
    E_Int jp = nj[i];
    E_Int kp = nk[i];
    readInt(s, si, 
            ptrFile, _convertEndian, si, si2);
    E_Int sizeBytes = 5*sizeof(E_Float)*ip*jp*kp;
    if (s != sizeBytes) // sep 5 * taille float * taille ndom
    {
      printf("Warning: binary file corrupted.\n");
      fclose(ptrFile);
      return 1;
    }
    
    E_Int start = 1;
    if (xyzdone == false)
    {
      FldArrayF* f = new FldArrayF(ip*jp*kp, 5);
      field.push_back(f);
    }
    else
    {
      FldArrayF* f = field[i];
      E_Int size = f->getSize();
      E_Int nfld = f->getNfld();
      f->reAlloc(size, 5+nfld);
      if (nfld == 4)
      {
        // Deplacement de cellN a la fin
        E_Float* pf = f->begin();
        for (E_Int p = 0; p < size; p++)
          pf[8*size+p] = pf[3*size+p];
      }
      start = nfld;
    }
    fread(field[i]->begin(start), sizeof(E_Float), 5*ip*jp*kp, ptrFile);
    fread(&ib, si, 1, ptrFile); // sep taille dom * 5 * taille float
  }
  fclose(ptrFile);

  // Build varString
  if (xyzdone == false)
    strcpy(varString, "ro,rou,rov,row,roE");
  else
  {
    if (field[0]->getNfld() == 9)
    {
      strcpy(varString, "x,y,z,ro,rou,rov,row,roE,cellN");
    }
    else
      strcpy(varString, "x,y,z,ro,rou,rov,row,roE");
  }
  return 0;
}

//=============================================================================
E_Int K_IO::GenIO::p3dfread(
  char* file, char*& varString,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& field)
{
  FILE* ptrFile;
  int s;
  int ib;
  vector<E_Int> nvar;

  // Un fichier a deja ete lu
  E_Boolean xyzdone = false;
  if (ni.size() != 0)
    xyzdone = true;

  // Le cellN a ete lu
  E_Boolean cellN = false;
  if (xyzdone == true)
  {
    E_Int l = strlen(varString);
    if (l >= 5 && varString[l-5] == 'c' && varString[l-4] == 'e' &&
        varString[l-3] == 'l' && varString[l-2] == 'l' &&
        varString[l-1] == 'N')
      cellN = true;
  }
  
  /* Constantes */
  E_Int si = sizeof(int);
  E_Int si2 = sizeof(int);
  
  // Check endianess
  E_Int ret = v3dCheckEndian(file); // similar to v3d
  if (ret == -1)
  {
    printf("Warning: fread: can not open file %s.\n", file);
    return 1;
  }
  if (ret == -2)
  {
    printf("Warning: fread: format in file %s is not recognised.\n", file);
    return 1;
  }
  else if (ret == 0)
  {
    _convertEndian = false;
    _intLength = 4;
  }
  else if (ret == 1)
  {
    _convertEndian = false;
    _intLength = 8;
  }
  else if (ret == 2)
  {
    _convertEndian = true;
    _intLength = 4;
  }
  else
  {
    _convertEndian = true;
    _intLength = 8;
  }

  /* Opening */
  ptrFile = fopen(file, "rb");

  // Read number of zones in file
  fread(&ib, si, 1, ptrFile); // sep
  E_Int nzones;
  if (_intLength == 4)
  {
    fread(&ib, si, 1, ptrFile);
    if (_convertEndian == true)
      nzones = IBE(ib);
    else
      nzones = ib;
  }
  else
  {
    fread(&s, si2, 1, ptrFile);
    if (_convertEndian == true)
      nzones = LBE(s);
    else
      nzones = s;
  }

  // Coherence
  if (xyzdone == true)
  {
    E_Int niSize = ni.size();
    if (niSize != nzones)
      printf("Warning: fread: the number of zones in f and xyz or q files is different.\n");
  }

  fread(&ib, si, 1, ptrFile); // sep taille int
  fread(&ib, si, 1, ptrFile); // sep nombre de zones * 3 * taille int
  
  // Read dimensions of zones
  for (E_Int i = 0; i < nzones; i++)
  {
    E_Int ip;
    if (_intLength == 4)
    {
      fread(&ib, si, 1, ptrFile);
      if (_convertEndian == true) ip = IBE(ib);
      else ip = ib;
    }
    else
    {
      fread(&s, si2, 1, ptrFile);
      if (_convertEndian == true) ip = LBE(s);
      else ip = s;
    } 
    
    E_Int jp;
    if (_intLength == 4)
    {
      fread(&ib, si, 1, ptrFile);
      if (_convertEndian == true) jp = IBE(ib);
      else jp = ib;
    }
    else
    {
      fread(&s, si2, 1, ptrFile);
      if (_convertEndian == true) jp = LBE(s);
      else jp = s;
    }
    
    E_Int kp;
    if (_intLength == 4)
    {
      fread(&ib, si, 1, ptrFile);
      if (_convertEndian == true) kp = IBE(ib);
      else kp = ib;
    }
    else
    {
      fread(&s, si2, 1, ptrFile);
      if (_convertEndian == true) kp = LBE(s);
      else kp = s;
    }

    E_Int nv;
    if (_intLength == 4)
    {
      fread(&ib, si, 1, ptrFile);
      if (_convertEndian == true) nv = IBE(ib);
      else nv = ib;
    }
    else
    {
      fread(&s, si2, 1, ptrFile);
      if (_convertEndian == true) nv = LBE(s);
      else nv = s;
    }
    nvar.push_back(nv);

    // Coherence
    if (xyzdone == true)
    {
      if (ip != ni[i] || jp != nj[i] || kp != nk[i])
        printf("Warning: zones in f file and xyz or q files are different.\n");
    }
    else
    {
      ni.push_back(ip); nj.push_back(jp); nk.push_back(kp);
    }
  }
  
  // Sep
  //E_Int sep;
  fread(&ib, si, 1, ptrFile); // sep nombre de zone * taille int * 4
  //if (_convertEndian == true) sep = IBE(ib);
  //else sep = ib;

  // Read additional variables
  for (E_Int i = 0; i < nzones; i++)
  {
    
    fread(&ib, si, 1, ptrFile); // sep taille dom * nvar * taille float
    //if (_convertEndian == true) sep = IBE(ib);
    //else sep = ib;

    E_Int ip = ni[i];
    E_Int jp = nj[i];
    E_Int kp = nk[i];
    E_Int nv = nvar[i];
    E_Int start = 1;
    if (xyzdone == false)
    {
      FldArrayF* f = new FldArrayF(ip*jp*kp, 5);
      field.push_back(f);
    }
    else
    {
      FldArrayF* f = field[i];
      E_Int size = f->getSize();
      E_Int nfld = f->getNfld();
      f->reAlloc(size, nv+nfld);
      if (cellN == true)
      {
        // Deplacement de cellN a la fin
        E_Float* pf = f->begin();
        for (E_Int p = 0; p < size; p++)
          pf[(nfld-1)*size+p] = pf[(nfld+nv-1)*size+p]; 
      }
      start = nfld;
    }
    fread(field[i]->begin(start), sizeof(E_Float), nv*ip*jp*kp, ptrFile);
    fread(&ib, si, 1, ptrFile); // sep taille dom * 3 * taille float
    //if (_convertEndian == true) sep = IBE(ib);
    //else sep = ib;
  }
  fclose(ptrFile);

  // Build varString
  char no[128];
  if (cellN == false)
  {
    E_Int nv = nvar[0];
    for (int p = 0 ; p < nv; p++)
    {
      if (strlen(varString) != 0) strcat(varString, ",");
      strcat(varString, "var");
      sprintf(no, "%d", p);
      strcat(varString, no);
    }
  }
  else
  {
    E_Int nv = nvar[0];
    E_Int l = strlen(varString);
    varString[l-6] = '\0';
    for (int p = 0 ; p < nv; p++)
    {
      if (strlen(varString) != 0)
        strcat(varString, ",");
      strcat(varString, "var");
      sprintf(no, "%d", p);
      strcat(varString, no);
    }
    strcat(varString, ",cellN");
  }

  return 0;
}

//=============================================================================
/* plot3dread.
   Read binary plot3d file as generated by ifort compiler. */
//=============================================================================
E_Int K_IO::GenIO::plot3dread(
  char* file, char*& varString,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& field, vector<char*>& zoneNames)
{
  varString = new char [K_ARRAY::VARSTRINGLENGTH];
  varString[0] = '\0';

  // Find the openable files
  char* xyzFile = NULL;
  char* qFile = NULL;
  char* fFile = NULL;
  E_Int res = findOpenablePlot3dFiles(file, xyzFile, qFile, fFile);
  
  if (res == 0)
  {
    return 1; // no openable file
  }

  // Read xyzFile if possible
  if (xyzFile != NULL)
  {
    res = xyzread(xyzFile, varString, ni, nj, nk,
                  field);
    delete [] xyzFile;
  }

  // Cree les noms de zones
  E_Int fieldSize = field.size();
  for (E_Int i = 0; i < fieldSize; i++)
  {
    char* zoneName = new char [128];
    sprintf(zoneName, "Zone%d", i);
    zoneNames.push_back(zoneName);
  }

  // Read q file if possible
  if (qFile != NULL)
  {
    res = p3dqread(qFile, varString, ni, nj, nk,
                   field);
    delete [] qFile;
  }

  // Read f file if possible
  if (fFile != NULL)
  {
    res = p3dfread(fFile, varString, ni, nj, nk,
                   field);
    delete [] fFile;
  }

  // Endian conversion of field
  if (_convertEndian == true && fieldSize > 0)
  {
    for (E_Int i = 0; i < fieldSize; i++) convertEndianField(*field[i]);
  } 
  return 0;
}

//=============================================================================
// xyzwrite
//=============================================================================
E_Int K_IO::GenIO::xyzwrite(
  char* file, char* varString, 
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& field, E_Int isize, E_Int rsize, 
  E_Boolean convertEndian)
{
  FILE* ptrFile;
  E_Int res, resCellN;

  /* Constants */
  E_Int si = sizeof(int);
  E_Int si2 = sizeof(E_LONG);

  // Open file
  ptrFile = fopen(file, "wb");
  if (ptrFile == NULL)
  {
    printf("Warning: xyzwrite: can not open file %s\n.", file);
    return 1;
  }

  // Get number of variables
  if (field.size() == 0) return 0; // nothing to write
  E_Int nvar = field[0]->getNfld();

  // Sep: taille int
  writeInt(isize, si,
           ptrFile, convertEndian, si, si2);

  // Nombre de zones
  E_Int nzones = field.size();
  writeInt(nzones, isize,
           ptrFile, convertEndian, si, si2);

  // Sep: taille int
  writeInt(isize, si,
           ptrFile, convertEndian, si, si2);

  // Sep: nombre de zones * 3 * taille int
  writeInt(nzones * isize * 3, si,
           ptrFile, convertEndian, si, si2);

  for (E_Int i = 0; i < nzones; i++)
  { 
    writeInt(ni[i], isize,
             ptrFile, convertEndian, si, si2);
    writeInt(nj[i], isize,
             ptrFile, convertEndian, si, si2);
    writeInt(nk[i], isize,
             ptrFile, convertEndian, si, si2);
  }
  
  // Sep : nombre de zones * 3 * taille int
  writeInt(nzones * isize * 3, si,
           ptrFile, convertEndian, si, si2);

  resCellN = isNamePresent("cellN", varString);
  if (resCellN == -1) resCellN = isNamePresent("cellnf", varString);
  if (resCellN == -1) resCellN = isNamePresent("cellNF", varString);

  nvar = 3 * sizeof(E_Float);
  if (resCellN != -1)
    nvar = nvar + isize;

  for (E_Int i = 0; i < nzones; i++)
  {
    // Get buffer
    FldArrayF buf(*field[i]);
    if (convertEndian == true)
      convertEndianField(buf);

    // Sep: taille dom * nvar * taille float
    writeInt(ni[i]*nj[i]*nk[i] * nvar, si,
             ptrFile, convertEndian, si, si2);

    res = isCoordinateXPresent(varString);
    fwrite(buf.begin(res+1), sizeof(E_Float), 
           ni[i]*nj[i]*nk[i], ptrFile);
    res = isCoordinateYPresent(varString);
    fwrite(buf.begin(res+1), sizeof(E_Float), 
           ni[i]*nj[i]*nk[i], ptrFile);
    res = isCoordinateZPresent(varString);
    fwrite(buf.begin(res+1), sizeof(E_Float), 
           ni[i]*nj[i]*nk[i], ptrFile);
    
    if (resCellN != -1)
    {
      if (isize == si)
      {
        int* iblank = new int[ni[i]*nj[i]*nk[i]];
        E_Float* pf = buf.begin(resCellN+1);
        for (E_Int ind = 0; ind < ni[i]*nj[i]*nk[i]; ind++)
          iblank[ind] = (int)pf[ind];
        fwrite(iblank, isize, ni[i]*nj[i]*nk[i], ptrFile);
        delete [] iblank;
      }
      else
      {
        E_LONG* iblank = new E_LONG[ni[i]*nj[i]*nk[i]];
        E_Float* pf = buf.begin(resCellN+1);
        for (E_Int ind = 0; ind < ni[i]*nj[i]*nk[i]; ind++)
          iblank[ind] = (E_LONG)pf[ind];
        fwrite(iblank, isize, ni[i]*nj[i]*nk[i], ptrFile);
        delete [] iblank;
      }
    }

    // Sep: taille dom * nvar * taille float
    writeInt(ni[i]*nj[i]*nk[i] * nvar, si,
             ptrFile, convertEndian, si, si2);
  }

  fclose(ptrFile);
  return 0;  
}

//=============================================================================
// qwrite
//=============================================================================
E_Int K_IO::GenIO::qwrite(
  char* file, char* varString, 
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& field, E_Int isize, E_Int rsize, 
  E_Boolean convertEndian)
{
  FILE* ptrFile;
  E_Int res;

  /* Constants */
  E_Int si = sizeof(int);
  E_Int si2 = sizeof(E_LONG);

  // Open file
  ptrFile = fopen(file, "wb");
  if (ptrFile == NULL)
  {
    printf("Warning: qwrite: can not open file %s\n.", file);
    return 1;
  }

  // Get number of variables
  if (field.size() == 0) return 0; // nothing to write
  E_Int nvar = field[0]->getNfld();

  // Sep: taille int
  writeInt(isize, si, ptrFile, convertEndian, si, si2);

  // Nombre de zones
  E_Int nzones = field.size();
  writeInt(nzones, isize,
           ptrFile, convertEndian, si, si2);

  // Sep: taille int
  writeInt(isize, si,
           ptrFile, convertEndian, si, si2);

  // Sep: nombre de zones * 3 * taille int
  writeInt(nzones * isize * 3, si,
           ptrFile, convertEndian, si, si2);

  for (E_Int i = 0; i < nzones; i++)
  { 
    writeInt(ni[i], isize,
             ptrFile, convertEndian, si, si2);
    writeInt(nj[i], isize,
             ptrFile, convertEndian, si, si2);
    writeInt(nk[i], isize,
             ptrFile, convertEndian, si, si2);
  }
  
  // Sep: nombre de zones * 3 * taille int
  writeInt(nzones*isize*3, si,
           ptrFile, convertEndian, si, si2);

  nvar = 5;

  for (E_Int i = 0; i < nzones; i++)
  {
    // Get buffer
    FldArrayF buf(*field[i]);
    if (convertEndian == true) convertEndianField(buf);
    
    // Sep : 4 * taille double
    writeInt(4 * sizeof(E_Float), si,
             ptrFile, convertEndian, si, si2);

    // Mach,...
    E_Float value = 0.;
    fwrite(&value, sizeof(E_Float), 1, ptrFile);
    fwrite(&value, sizeof(E_Float), 1, ptrFile);
    fwrite(&value, sizeof(E_Float), 1, ptrFile);
    fwrite(&value, sizeof(E_Float), 1, ptrFile);

    // Sep : 4 * taille double
    writeInt(4 * sizeof(E_Float), si,
             ptrFile, convertEndian, si, si2);

    // Sep : taille dom * 5 * taille float
    writeInt(ni[i]*nj[i]*nk[i] * sizeof(E_Float) * nvar, si,
             ptrFile, convertEndian, si, si2);

    res = isDensityPresent(varString);
    fwrite(buf.begin(res+1), sizeof(E_Float), 
           ni[i]*nj[i]*nk[i], ptrFile);
    res = isMomentumXPresent(varString);
    fwrite(buf.begin(res+1), sizeof(E_Float), 
           ni[i]*nj[i]*nk[i], ptrFile);
    res = isMomentumYPresent(varString);
    fwrite(buf.begin(res+1), sizeof(E_Float), 
           ni[i]*nj[i]*nk[i], ptrFile);
    res = isMomentumZPresent(varString);
    fwrite(buf.begin(res+1), sizeof(E_Float), 
           ni[i]*nj[i]*nk[i], ptrFile); 
    res = isEnergyStagnationDensityPresent(varString);
    fwrite(buf.begin(res+1), sizeof(E_Float), 
           ni[i]*nj[i]*nk[i], ptrFile);

    // Sep : taille dom * nvar * taille float
    writeInt(ni[i]*nj[i]*nk[i] * sizeof(E_Float) * nvar, si,
             ptrFile, convertEndian, si, si2);
  }

  fclose(ptrFile);
  return 0;
}

//=============================================================================
// fwrite
// Ecrit les autres champs (autre que x,y,z)
// Si writeConservative=True, ecrit aussi les variables conservatives
//=============================================================================
E_Int K_IO::GenIO::ffwrite(
  char* file, char* varString, 
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& field, E_Int isize, E_Int rsize, 
  E_Boolean convertEndian, E_Boolean writeConservative)
{
  FILE* ptrFile;

  /* Constants */
  E_Int si = sizeof(int);
  E_Int si2 = sizeof(E_LONG);

  // Get number of variables
  if (field.size() == 0) return 0; // nothing to write
  E_Int nvar = field[0]->getNfld();

  // Open file
  ptrFile = fopen(file, "wb");
  if (ptrFile == NULL)
  {
    printf("Warning: ffwrite: can not open file %s\n.", file);
    return 1;
  }

  // Sep : taille int
  writeInt(isize, si,
           ptrFile, convertEndian, si, si2);

  // Nombre de zones
  E_Int nzones = field.size();
  writeInt(nzones, isize,
           ptrFile, convertEndian, si, si2);

  // Sep : taille int
  writeInt(isize, si,
           ptrFile, convertEndian, si, si2);

  // Sep : nombre de zones * 3 * taille int
  writeInt(nzones * isize * 3, si,
           ptrFile, convertEndian, si, si2);

  for (E_Int i = 0; i < nzones; i++)
  { 
    writeInt(ni[i], isize,
             ptrFile, convertEndian, si, si2);
    writeInt(nj[i], isize,
             ptrFile, convertEndian, si, si2);
    writeInt(nk[i], isize,
             ptrFile, convertEndian, si, si2);
  }
  
  // Sep : nombre de zones * 3 * taille int
  writeInt(nzones * isize * 3, si,
           ptrFile, convertEndian, si, si2);

  for (E_Int i = 0; i < nzones; i++)
  {
    nvar = field[i]->getNfld();

    // Get buffer
    FldArrayF buf(*field[i]);
    if (convertEndian == true) convertEndianField(buf);
    
    // Sep : 4 * taille double
    writeInt(4 * sizeof(E_Float), si,
             ptrFile, convertEndian, si, si2);

    // Mach,...
    E_Float value = 0.;
    fwrite(&value, sizeof(E_Float), 1, ptrFile);
    fwrite(&value, sizeof(E_Float), 1, ptrFile);
    fwrite(&value, sizeof(E_Float), 1, ptrFile);
    fwrite(&value, sizeof(E_Float), 1, ptrFile);

    // Sep : 4 * taille double
    writeInt(4 * sizeof(E_Float), si,
             ptrFile, convertEndian, si, si2);

    // Sep : taille dom * 5 * taille float
    writeInt(ni[i]*nj[i]*nk[i] * sizeof(E_Float) * nvar, si,
             ptrFile, convertEndian, si, si2);

    // Position des variables 
    E_Int resx = isCoordinateXPresent(varString);
    E_Int resy = isCoordinateYPresent(varString);
    E_Int resz = isCoordinateZPresent(varString);
    E_Int res1 = isDensityPresent(varString);
    E_Int res2 = isMomentumXPresent(varString);
    E_Int res3 = isMomentumYPresent(varString);
    E_Int res4 = isMomentumZPresent(varString);
    E_Int res5 = isEnergyStagnationDensityPresent(varString);

    for (E_Int nv = 1; nv <= nvar; nv++)
    {
      if (writeConservative == false)
      {
        if (nv != resx && nv != resy && nv != resz &&
            nv != res1 && nv != res2 && nv != res3 &&
            nv != res4 && nv != res5)
          fwrite(buf.begin(nv), sizeof(E_Float), 
                 ni[i]*nj[i]*nk[i], ptrFile);
      }
      else
      {
        if (nv != resx && nv != resy && nv != resz)
          fwrite(buf.begin(nv), sizeof(E_Float), 
                 ni[i]*nj[i]*nk[i], ptrFile);
      }
    }

    // Sep : taille dom * nvar * taille float
    writeInt(ni[i]*nj[i]*nk[i] * sizeof(E_Float) * nvar, si,
             ptrFile, convertEndian, si, si2);
  }

  fclose(ptrFile);
  return 0;  
}

//=============================================================================
// plot3dwrite
// write plot3d file, with options
// Return 1 upon failure.
//=============================================================================
E_Int K_IO::GenIO::plot3dwrite(
  char* file, char* dataFmt, char* varString, 
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& field,
  vector<char*>& zoneNames, E_Int isize, E_Int rsize, 
  E_Boolean convertEndian)
{
  // Find writable files
  char* xyzFile = NULL;
  char* qFile = NULL;
  char* fFile = NULL;
  findWritablePlot3dFiles(file, varString,
                          xyzFile, qFile, fFile);

  if (xyzFile != NULL)
  {
    E_Int resx = isCoordinateXPresent(varString);
    E_Int resy = isCoordinateYPresent(varString);
    E_Int resz = isCoordinateZPresent(varString);
    if (resx == -1 || resy == -1 || resz == -1)
    {
      if (resx == -1)
        printf("Warning: xyz can not be written because no x variable was found in array.\n");
      if (resy == -1)
        printf("Warning: xyz can not be written because no y variable was found in array.\n");
      if (resz == -1)
        printf("Warning: xyz can not be written because no z variable was found in array.\n");
    }
    else
    {
      xyzwrite(xyzFile, varString,
               ni, nj, nk, field, isize, rsize, convertEndian);
    }
    delete [] xyzFile;
  }

  if (qFile != NULL)
  {
    E_Int resro = isDensityPresent(varString);
    E_Int resrou = isMomentumXPresent(varString);
    E_Int resrov = isMomentumYPresent(varString);
    E_Int resrow = isMomentumZPresent(varString);
    E_Int resroE = isEnergyStagnationDensityPresent(varString);
    if (resro == -1 || resrou == -1 || resrov == -1 
        || resrow == -1 || resroE == -1)
    {
      if (resro == -1)
        printf("Warning: q file can not be written because no ro variable was found in array.\n");
      if (resrou == -1)
        printf("Warning: q file can not be written because no rou variable was found in array.\n");
      if (resrov == -1)
        printf("Warning: q file can not be written because no rov variable was found in array.\n");
      if (resrow == -1)
        printf("Warning: q file can not be written because no row variable was found in array.\n");
      if (resroE == -1)
        printf("Warning: q file can not be written because no roE variable was found in array.\n");
    }
    else
    {
      qwrite(qFile, varString,
             ni, nj, nk, field, isize, rsize, convertEndian); 
    }
    delete [] qFile;
  }
 
  if (fFile != NULL)
  {
    E_Boolean writeConservative = true;
    if (qFile == NULL) writeConservative = false;
    ffwrite(fFile, varString,
            ni, nj, nk, field, isize, rsize, convertEndian, 
            writeConservative); 
    delete [] fFile;
  }

  return 0;

}
