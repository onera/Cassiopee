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

// Binary archive (CEDRE) file support

#define SENTINELLE -1.79769e+308
#define STRINGLENGTH 256

# include "GenIO.h"
#include "kcore.h"
# include <stdio.h>
# include <string.h>
# include "Array/Array.h"
# include <map>

// for dbx
#include <iostream>

using namespace std;
using namespace K_FLD;

class mesh
{

public:
  mesh() { _isomu=NULL; _ifacu=NULL; _icelu=NULL; _pe=NULL; _nbpoints=NULL; 
    _numpoints=NULL; _x=NULL; _y=NULL; _z=NULL; };
  ~mesh() { delete [] _isomu; delete [] _ifacu; delete [] _icelu; 
    delete [] _pe; delete [] _nbpoints; delete [] _numpoints; 
    delete [] _x; delete [] _y; delete [] _z;
    for (size_t i = 0; i < _nfields.size(); i++) delete [] _nfields[i];
    for (size_t i = 0; i < _cfields.size(); i++) delete [] _cfields[i];
    for (size_t i = 0; i < _nfieldNames.size(); i++) delete [] _nfieldNames[i];
    for (size_t i = 0; i < _cfieldNames.size(); i++) delete [] _cfieldNames[i];
  };
  
public:
  // IN: Numero absolu + temps permettent de reperer un bloc de facon unique
  // numero absolu du sous domaine
  unsigned int _nd;
  // temps
  double _temps;
  
  // numero du domaine utilisateur auquel ce bloc appartient
  unsigned int _numutil;
  // nom du du domaine
  char _nom[STRINGLENGTH];
  // type du maillage
  unsigned char _type;
  // situation des variables (0,1)
  unsigned char _situ;
  // numero du melange
  int _num_mel;
  // numero du groupe
  int _num_grp;
  // nombre de sommets
  int _nsom;
  // nombre de faces
  int _nfac;
  // nombre de cellules limites
  int _nflm;
  int _nflp;
  // nombre de cellules internes
  int _nceli;
  // isomu
  int* _isomu;
  // ifacu
  int* _ifacu;
  // icelu
  int* _icelu;
  // elmin
  unsigned int _elminx;
  // nelem
  unsigned int _nelemx;
  // parent elements
  unsigned int* _pe;
  // nbpoints - nbre de points par face
  unsigned char* _nbpoints;
  // numpoints - connectivite face->noeuds
  unsigned int* _numpoints;
  // imax jmax kmax
  unsigned short _imax;
  unsigned short _jmax;
  unsigned short _kmax;
  // elmin
  unsigned short _elmin;
  // nelem
  unsigned short _nelem;
  // coords
  double* _x;
  double* _y;
  double* _z;
  // node fields
  vector<double*> _nfields;
  // node fields names
  vector<char*> _nfieldNames;
  // center fields
  vector<double*> _cfields;
  // center fields names
  vector<char*> _cfieldNames;
  };

// Verifie si le domaine nd existe dans meshes, sinon le cree
// il faudrait aussi pouvoir utiliser le temps
mesh* checkMesh(unsigned int nd, std::map<unsigned int,mesh*>& meshes)
{
  std::map<unsigned int,mesh*>::iterator it;
  it = meshes.find(nd);
  
  if (it != meshes.end()) return it->second;
  else
  {
    mesh* m = new mesh();
    meshes[nd] = m;
    return m;
  }
}

// Read a block name jusqu'au caractere nul
void readBlockName(FILE* ptrFile, char* name)
{
  int c;
  E_Int i = 0;
  while ((c = fgetc(ptrFile)) != '\0')
  { name[i] = c; i++; }
  name[i] = '\0';
  printf("blockname: %s\n", name);
}

void readNElemElmin(FILE* ptrFile, unsigned int& nelem, unsigned int& elmin)
{
  fread(&elmin, sizeof(unsigned int), 1, ptrFile); elmin = UIBE(elmin);
  printf("elmin=%u\n", elmin);
  fread(&nelem, sizeof(unsigned int), 1, ptrFile); nelem = UIBE(nelem);
  printf("nelem=%u\n", nelem);
}

// Lit et decompresse un champ de floats
double* readAndUncompress(FILE* ptrFile, unsigned int nelem, unsigned int elmin)
{
  // Lit comp
  unsigned char comp;
  fread(&comp, sizeof(unsigned char), 1, ptrFile);
  printf("compression=%hhu\n", comp);
  double vmin = 0.;
  if (comp != 0) { fread(&vmin, sizeof(double), 1, ptrFile); vmin = DBE(vmin); }
  double vmax = 0.;
  if (comp != 0 && comp != 4) { fread(&vmax, sizeof(double), 1, ptrFile); vmax = DBE(vmax); }
  
  if (comp == 0) // pas de compression
  {
    double* F = new double [nelem];
    fread(F, sizeof(double), nelem, ptrFile);
    for (size_t i = 0; i < nelem; i++) F[i] = DBE(F[i]);
    return F;
  }
  else if (comp == 1)
  {
    unsigned int* F = new unsigned int [nelem];
    fread(F, sizeof(unsigned int), nelem, ptrFile);
    for (size_t i = 0; i < nelem; i++) F[i] = UIBE(F[i]);
    double* F2 = new double [nelem];
    //double fmin = 0;
    double fmax = 4294967295.;
    for (size_t i = 0; i < nelem; i++) F2[i] = vmin + F[i]*(vmax-vmin)/fmax;
    delete [] F;
    return F2;
  }
  else if (comp == 2)
  {
    unsigned short* F = new unsigned short [nelem];
    fread(F, sizeof(unsigned short), nelem, ptrFile); 
    double* F2 = new double [nelem];
    //double fmin = 0;
    double fmax = 65535.;
    for (size_t i = 0; i < nelem; i++) F2[i] = vmin + F[i]*(vmax-vmin)/fmax;
    delete [] F;
    return F2;
  }
  else if (comp == 3)
  {
    unsigned char* F = new unsigned char [nelem];
    fread(F, sizeof(unsigned char), nelem, ptrFile); 
    double* F2 = new double [nelem];
    //double fmin = 0;
    double fmax = 255.;
    for (size_t i = 0; i < nelem; i++) F2[i] = vmin + F[i]*(vmax-vmin)/fmax;
    delete [] F; 
    return F2;
  }
  else if (comp == 4)
  {
    double* F = new double [nelem];
    for (size_t i = 0; i < nelem; i++) F[i] = vmin;
    return F;
  }
  else if (comp == 5)
  {
    float* F = new float [nelem];
    fread(F, sizeof(float), nelem, ptrFile); 
    double* F2 = new double [nelem];
    //double fmin = 0;
    double fmax = K_CONST::E_MAX_FLOAT;
    for (size_t i = 0; i < nelem; i++) F2[i] = vmin + F[i]*(vmax-vmin)/fmax;
    delete [] F;
    return F2;
  }
  return NULL;
}

// Lit la version du fichier
void readVersion(FILE* ptrFile, unsigned char* version)
{
  int c;
  E_Int i = 0;
  while ((c = fgetc(ptrFile)) != '\0') { version[i] = c; i++; }
  version[i] = '\0';
  printf("version: %u\n", *version);
}

// Lit le titre du fichier
void readTitle(FILE* ptrFile, char* titre, char* date, char* machine, 
               double& tempsZero, double& epsilon)
{
  int c; E_Int i;
  if (titre != NULL)
  {
    i = 0;
    while ((c = fgetc(ptrFile)) != '\0' && i < STRINGLENGTH-1) { titre[i] = c; i++; }
    titre[i] = '\0';
  }
  i = 0;
  while ((c = fgetc(ptrFile)) != '\0' && i < STRINGLENGTH-1) { date[i] = c; i++; }
  date[i] = '\0';
  i = 0;
  while ((c = fgetc(ptrFile)) != '\0' && i < STRINGLENGTH-1) { machine[i] = c; i++; }
  machine[i] = '\0';
  fread(&tempsZero, sizeof(double), 1, ptrFile); tempsZero = DBE(tempsZero);
  fread(&epsilon, sizeof(double), 1, ptrFile); epsilon = DBE(epsilon);
  
  if (titre != NULL) printf("titre: %s\n", titre);
  printf("date: %s\n", date);
  printf("machine: %s\n", machine);
  printf("temps0: %g\n", tempsZero);
  printf("epsilon: %g\n", epsilon);
}

// Lit les donnes globales
void readGlobal(FILE* ptrFile, unsigned char& solverType, unsigned char& dimField)
{
  fread(&solverType, sizeof(unsigned char), 1, ptrFile);
  fread(&dimField, sizeof(unsigned char), 1, ptrFile);
  printf("solverType: %u\n", solverType);
  printf("dimField: %u\n", dimField);
}

// Lit une espece
void readEspece(FILE* ptrFile, int& numMel, char* nomMel, int& nespeces, char* nomEspece)
{
  int c;
  // numero du melange
  fread(&numMel, sizeof(int), 1, ptrFile); numMel = IBE(numMel);
  E_Int i = 0;
  while ((c = fgetc(ptrFile) && i < STRINGLENGTH-1) != '\0') { nomMel[i] = c; i++; }
  nomMel[i] = '\0';
  fread(&nespeces, sizeof(int), 1, ptrFile); nespeces = IBE(nespeces);
  i = 0;
  for (E_Int j = 0; j < nespeces; j++)
  {
    while ((c = fgetc(ptrFile) && i < STRINGLENGTH-1) != '\0') { nomEspece[i] = c; i++; }
    nomEspece[i] = '\0';
  }
  printf("numMel: %d\n", numMel);
  printf("nomMel: %s\n", nomMel);
  printf("nbreEspeces: %d\n", nespeces);
  for (E_Int i = 0; i < nespeces; i++) printf(SF_D_ " %s\n", i, nomEspece);
}

// Lit les donnees scalaires
void readScalar(FILE* ptrFile, int& numGrp, char* nomGrp, int& nelem, 
                char* nomScalar, unsigned char& typeSca)
{
  int c;
  fread(&numGrp, sizeof(int), 1, ptrFile); numGrp = IBE(numGrp);
  E_Int i = 0;
  while ((c = fgetc(ptrFile)) != '\0') { nomGrp[i] = c; i++; }
  nomGrp[i] = '\0';
  fread(&nelem, sizeof(int), 1, ptrFile); nelem = IBE(nelem);
  i = 0;
  while ((c = fgetc(ptrFile)) != '\0' && i < STRINGLENGTH-1) { nomScalar[i] = c; i++; }
  nomScalar[i] = '\0';
  fread(&typeSca, sizeof(unsigned char), 1, ptrFile);
}

// Lit une unite
void readUnit(FILE* ptrFile, char* name, char* unit)
{
  int c;
  E_Int i = 0;
  while ((c = fgetc(ptrFile)) != '\0' && i < STRINGLENGTH-1) { name[i] = c; i++; }
  name[i] = '\0';
  printf("name %s\n", name);
  i = 0;
  while ((c = fgetc(ptrFile)) != '\0' && i < STRINGLENGTH-1) { unit[i] = c; i++; }
  unit[i] = '\0';
  printf("unit %s\n", unit);
}

// Lit les donnees thermo
void readThermo(FILE* ptrFile, unsigned int& read, char* name,  
  unsigned int& nGe, unsigned int& nGr,
  unsigned int& ne, unsigned int& nr, 
  unsigned int*& dthGe, unsigned int*& dthGr,
  unsigned int*& dthe, double*& dthr)
{
  fread(&read, sizeof(unsigned int), 1, ptrFile); read = UIBE(read);
  printf("thermo: read %u\n", read);
  
  int c;
  E_Int i = 0;
  while ( (c = fgetc(ptrFile)) != '\0' && i < STRINGLENGTH-1) { name[i] = c; i++; }
  name[i] = '\0';
  printf("thermo name : %s\n", name);
  if (read == 1) return;
  fread(&nGe, sizeof(unsigned int), 1, ptrFile); nGe = UIBE(nGe);
  fread(&nGr, sizeof(unsigned int), 1, ptrFile); nGr = UIBE(nGr);
  fread(&ne, sizeof(unsigned int), 1, ptrFile); ne = UIBE(ne);
  fread(&nr, sizeof(unsigned int), 1, ptrFile); nr = UIBE(nr);
  printf("nGe: %u %u %u %u\n", nGe, nGr, ne, nr);
  dthGe = new unsigned int [nGe];
  dthGr = new unsigned int [nGr];
  dthe = new unsigned int [ne];
  dthr = new double [nr];
  fread(&dthGe, sizeof(unsigned int), nGe, ptrFile); 
  fread(&dthGr, sizeof(unsigned int), nGr, ptrFile);
  fread(&dthe, sizeof(unsigned int), ne, ptrFile);
  fread(&dthr, sizeof(double), nr, ptrFile);
}
 
// Lit la structure d'un maillage
void readStructure(FILE* ptrFile, std::map<unsigned int, mesh*>& meshes)
{  
  unsigned int nd;
  fread(&nd, sizeof(unsigned int), 1, ptrFile); nd = UIBE(nd);
  mesh* m = checkMesh(nd, meshes);
  
  fread(&m->_numutil, sizeof(unsigned int), 1, ptrFile); m->_numutil = UIBE(m->_numutil);

  E_Int i = 0; int c;
  while ((c = fgetc(ptrFile)) != '\0' && i < STRINGLENGTH-1) { m->_nom[i] = c; i++; }
  m->_nom[i] = '\0';
  // type=0 (structure), type=1 (mixte), type=2 (elements finis), type=3 (particules)
  fread(&m->_type, sizeof(unsigned char), 1, ptrFile); 
  fread(&m->_situ, sizeof(unsigned char), 1, ptrFile); 
  fread(&m->_num_mel, sizeof(int), 1, ptrFile); m->_num_mel = IBE(m->_num_mel);
  fread(&m->_num_grp, sizeof(int), 1, ptrFile); m->_num_grp = IBE(m->_num_grp);
  printf("Structure: nd=%u numutil=%u\n", nd, m->_numutil);
  printf("nom dom=%s\n", m->_nom);
  printf("type=%u, situ=%u\n", m->_type, m->_situ);
}

void readDomutil(FILE* ptrFile, unsigned int& numUti, char* nomUti)
{
  int c;
  E_Int i = 0;
  fread(&numUti, sizeof(unsigned int), 1, ptrFile); numUti = UIBE(numUti);
  while ((c = fgetc(ptrFile)) != '\0' && i < STRINGLENGTH-1) { nomUti[i] = c; i++; }
  nomUti[i] = '\0';
  printf("nom du domaine utilisateur %s\n", nomUti);
}

void readNumerotation(FILE* ptrFile, std::map<unsigned int, mesh*>& meshes)
{
  unsigned int nd;
  fread(&nd, sizeof(unsigned int), 1, ptrFile); nd = UIBE(nd);
  mesh* m = checkMesh(nd, meshes);
  
  fread(&m->_temps, sizeof(double), 1, ptrFile); m->_temps = DBE(m->_temps);
  
  fread(&m->_nsom, sizeof(int), 1, ptrFile); m->_nsom = IBE(m->_nsom);
  fread(&m->_nfac, sizeof(int), 1, ptrFile); m->_nfac = IBE(m->_nfac);
  fread(&m->_nflm, sizeof(int), 1, ptrFile); m->_nflm = IBE(m->_nflm);
  fread(&m->_nflp, sizeof(int), 1, ptrFile); m->_nflp = IBE(m->_nflp);
  fread(&m->_nceli, sizeof(int), 1, ptrFile); m->_nceli = IBE(m->_nceli);
  printf("npts=%d nfaces=%d \n", m->_nsom, m->_nfac);
  printf("ncellLim=%d, ncellPart=%d, ncelli=%d\n",
         m->_nflm, m->_nflp, m->_nceli);
  m->_isomu = new int [m->_nsom];
  fread(m->_isomu, sizeof(int), m->_nsom, ptrFile);
  int* pt = m->_isomu;
  for (E_Int i = 0; i < m->_nsom; i++) pt[i] = IBE(pt[i]);

  m->_ifacu = new int [m->_nfac];
  fread(m->_ifacu, sizeof(int), m->_nfac, ptrFile);
  pt = m->_ifacu;
  for (E_Int i = 0; i < m->_nfac; i++) pt[i] = IBE(pt[i]);
  
  E_Int size = m->_nceli+m->_nflm+m->_nflp;
  m->_icelu = new int [size];
  fread(m->_icelu, sizeof(int), size, ptrFile);
  pt = m->_icelu;
  for (E_Int i = 0; i < size; i++) pt[i] = IBE(pt[i]);
}

void readVariables(FILE* ptrFile, unsigned int& numUti, unsigned int& nvar, char* nomVar,
  char* nomGrp, char* nomCat, unsigned char& glob, unsigned int& elmin,
  unsigned int& nelem, double*& data)
{
  fread(&numUti, sizeof(unsigned int), 1, ptrFile); numUti = UIBE(numUti);
  fread(&nvar, sizeof(unsigned int), 1, ptrFile); nvar = UIBE(nvar);
  
}

// lecture de la connectivite
void readConnexion(FILE* ptrFile, std::map<unsigned int, mesh*>& meshes)
{
  unsigned int nd;
  fread(&nd, sizeof(unsigned int), 1, ptrFile); nd = UIBE(nd);
  mesh* m = checkMesh(nd, meshes);
  
  fread(&m->_temps, sizeof(double), 1, ptrFile); m->_temps = DBE(m->_temps);
  fread(&m->_elminx, sizeof(unsigned int), 1, ptrFile); m->_elminx = UIBE(m->_elminx);
  fread(&m->_nelemx, sizeof(unsigned int), 1, ptrFile); m->_nelemx = UIBE(m->_nelemx);
  
  E_Int size = 2*m->_nelemx;
  m->_pe = new unsigned int [size];
  fread(m->_pe, sizeof(unsigned int), size, ptrFile);
  for (E_Int i = 0; i < size; i++) m->_pe[i] = UIBE(m->_pe[i]);
  
  m->_nbpoints = new unsigned char [m->_nelemx];
  fread(m->_nbpoints, sizeof(unsigned char), m->_nelemx, ptrFile);
  
  E_Int nbPointsTot = 0;
  for (size_t i = 0; i < m->_nelemx; i++) nbPointsTot += (E_Int)m->_nbpoints[i];
    
  size = nbPointsTot;
  m->_numpoints = new unsigned int [size];
  fread(m->_numpoints, sizeof(unsigned int), size, ptrFile);
  for (E_Int i = 0; i < size; i++) m->_numpoints[i] = UIBE(m->_numpoints[i]);
}

// Lit les coords maillage structure (type=0 ou type=1)
void readMaillage0(FILE* ptrFile, mesh* m)
{
  double temps;
  fread(&temps, sizeof(double), 1, ptrFile); temps = DBE(temps);
  
  fread(&m->_imax, sizeof(unsigned short), 1, ptrFile); 
  fread(&m->_jmax, sizeof(unsigned short), 1, ptrFile); 
  fread(&m->_kmax, sizeof(unsigned short), 1, ptrFile);
  
  m->_x = readAndUncompress(ptrFile, m->_imax*m->_jmax*m->_kmax, 0);
  m->_y = readAndUncompress(ptrFile, m->_imax*m->_jmax*m->_kmax, 0);
  m->_z = readAndUncompress(ptrFile, m->_imax*m->_jmax*m->_kmax, 0);
}

// Lit un maillage non structure (type=2 ou type=3)
// Lit les coordonnees
void readMaillage1(FILE* ptrFile, mesh* m)
{
  double temps;
  fread(&temps, sizeof(double), 1, ptrFile); temps = DBE(temps);
  printf("temps=" SF_F_ "\n", temps);
  
  unsigned int elmin; unsigned int nelem;
  readNElemElmin(ptrFile, nelem, elmin);
  
  m->_x = readAndUncompress(ptrFile, nelem, elmin);
  m->_y = readAndUncompress(ptrFile, nelem, elmin);
  m->_z = readAndUncompress(ptrFile, nelem, elmin);
}

// Lecture generale de maillage
void readMaillage(FILE* ptrFile, std::map<unsigned int, mesh*>& meshes)
{
  unsigned int nd;
  fread(&nd, sizeof(unsigned int), 1, ptrFile); nd = UIBE(nd);
  mesh* m = checkMesh(nd, meshes);
  printf("reading type=%d for mesh no %d\n", m->_type, nd);
  if (m->_type == 0 || m->_type == 1) readMaillage0(ptrFile, m);
  else readMaillage1(ptrFile, m);
}

// Lecture de CELLULES
// contient un champ en centres 0: cellule interne, 1: cellule frontiere
void readCellules(FILE* ptrFile, std::map<unsigned int, mesh*>& meshes)
{
  unsigned int nd;
  fread(&nd, sizeof(unsigned int), 1, ptrFile); nd = UIBE(nd);
  double temps;
  fread(&temps, sizeof(double), 1, ptrFile); temps = DBE(temps);
  //mesh* m = checkMesh(nd, meshes);
  
  unsigned int elmin; unsigned int nelem;
  readNElemElmin(ptrFile, nelem, elmin);
  
  unsigned char* type = new unsigned char [nelem];
  fread(type, sizeof(unsigned char), nelem, ptrFile);
  delete [] type;
}

// Lecture des centres de gravite des cellules
void readGraviteCellule(FILE* ptrFile, std::map<unsigned int, mesh*>& meshes)
{
  unsigned int nd;
  fread(&nd, sizeof(unsigned int), 1, ptrFile); nd = UIBE(nd);
  double temps;
  fread(&temps, sizeof(double), 1, ptrFile); temps = DBE(temps);
  //mesh* m = checkMesh(nd, meshes);
  
  unsigned int elmin; unsigned int nelem;
  readNElemElmin(ptrFile, nelem, elmin);
  
  double* X = readAndUncompress(ptrFile, nelem, elmin);
  delete [] X;
  double* Y = readAndUncompress(ptrFile, nelem, elmin);
  delete [] Y;
  double* Z = readAndUncompress(ptrFile, nelem, elmin);
  delete [] Z;
}

// Lecture des centres de gravite des cellules
void readGraviteFace(FILE* ptrFile, std::map<unsigned int, mesh*>& meshes)
{
  unsigned int nd;
  fread(&nd, sizeof(unsigned int), 1, ptrFile); nd = UIBE(nd);
  double temps;
  fread(&temps, sizeof(double), 1, ptrFile); temps = DBE(temps);
  //mesh* m = checkMesh(nd, meshes);
  
  unsigned int elmin; unsigned int nelem;
  readNElemElmin(ptrFile, nelem, elmin);
  
  double* X = readAndUncompress(ptrFile, nelem, elmin);
  delete [] X;
  double* Y = readAndUncompress(ptrFile, nelem, elmin);
  delete [] Y;
  double* Z = readAndUncompress(ptrFile, nelem, elmin);
  delete [] Z;
}

// Surface pour visu
void readSurface(FILE* ptrFile, std::map<unsigned int, mesh*>& meshes)
{
  unsigned int nd;
  fread(&nd, sizeof(unsigned int), 1, ptrFile); nd = UIBE(nd);
  double temps;
  fread(&temps, sizeof(double), 1, ptrFile); temps = DBE(temps);
  //mesh* m = checkMesh(nd, meshes);
  
  // numero absolu de la surface
  unsigned int no;
  fread(&no, sizeof(unsigned int), 1, ptrFile); no = UIBE(no);
  // nom de la surface
  char nomSurf[STRINGLENGTH];
  E_Int i = 0; char c;
  while ((c = fgetc(ptrFile)) != '\0' && i < STRINGLENGTH-1) { nomSurf[i] = c; i++; }
  nomSurf[i] = '\0';
  printf("nom surf=%s\n", nomSurf);
  // type de la surface
  unsigned char typeSurf;
  fread(&typeSurf, sizeof(unsigned char), 1, ptrFile);
  
  unsigned int elmin; unsigned int nelem;
  readNElemElmin(ptrFile, nelem, elmin);
  
  // no des faces de cette surface
  unsigned int* connect = new unsigned int [nelem];
  fread(connect, sizeof(unsigned int), nelem, ptrFile);
  delete [] connect;
  // type de ces faces (tag)
  unsigned char* type = new unsigned char [nelem];
  fread(type, sizeof(unsigned char), nelem, ptrFile);
  delete [] type; 
}

// Lecture d'un champ volumique
void readValeurVolumique(FILE* ptrFile, std::map<unsigned int, mesh*>& meshes)
{
  // champ moyen ou instantane
  unsigned char moyen;
  fread(&moyen, sizeof(unsigned char), 1, ptrFile); //moyen = UIBE(moyen);
  printf("Champ moyen=%hhu\n", moyen);
  
  unsigned int nd;
  fread(&nd, sizeof(unsigned int), 1, ptrFile); nd = UIBE(nd);
  double temps;
  fread(&temps, sizeof(double), 1, ptrFile); temps = DBE(temps);
  printf("temps=" SF_F_ "\n", temps);
  mesh* m = checkMesh(nd, meshes);

  unsigned char classe;
  fread(&classe, sizeof(unsigned char), 1, ptrFile); //classe = UIBE(classe);
  printf("Champ de classe=%hhu\n", classe);
  
  // Unknown thing that is here (not in spec)
  char toto[20];
  fread(toto, sizeof(char), 8, ptrFile);
  
  // nom du champ
  char* nomField = new char [STRINGLENGTH];
  E_Int i = 0; char c;
  while ((c = fgetc(ptrFile)) != '\0' && i < STRINGLENGTH-1) { nomField[i] = c; i++; }
  nomField[i] = '\0';
  
  unsigned int elmin; unsigned int nelem;
  readNElemElmin(ptrFile, nelem, elmin);
  
  printf("nom champ=%s\n", nomField);
  printf("taille %u (ncell=%d)\n", nelem, m->_nceli);
  
  double* F = readAndUncompress(ptrFile, nelem, elmin);
  
  if ((int)nelem == m->_nsom) { m->_nfieldNames.push_back(nomField); m->_nfields.push_back(F); }
  else if ((int)nelem >= m->_nceli) { m->_cfieldNames.push_back(nomField); m->_cfields.push_back(F); }
  else delete [] F;
}

// Lecture d'un vecteur volumique
void readVecteurVolumique(FILE* ptrFile, std::map<unsigned int, mesh*>& meshes)
{
  // champ moyen ou instantane
  unsigned char moyen;
  fread(&moyen, sizeof(unsigned char), 1, ptrFile); //moyen = UIBE(moyen);
  printf("Champ moyen=%hhu\n", moyen);
  
  unsigned int nd;
  fread(&nd, sizeof(unsigned int), 1, ptrFile); nd = UIBE(nd);
  double temps;
  fread(&temps, sizeof(double), 1, ptrFile); temps = DBE(temps);
  printf("temps=" SF_F_ "\n", temps);
  mesh* m = checkMesh(nd, meshes);

  unsigned char classe;
  fread(&classe, sizeof(unsigned char), 1, ptrFile); //classe = UIBE(classe);
  printf("Champ de classe=%hhu\n", classe);
  
  // Unknown thing that is here (not in spec)
  char toto[20];
  fread(toto, sizeof(char), 8, ptrFile);
  
  // nom du champ
  char nomField[STRINGLENGTH];
  E_Int i = 0; char c;
  while ((c = fgetc(ptrFile)) != '\0' && i < STRINGLENGTH-1) { nomField[i] = c; i++; }
  nomField[i] = '\0';
  unsigned int elmin; unsigned int nelem;
  readNElemElmin(ptrFile, nelem, elmin);
  
  printf("nom champ=%s\n", nomField);
  printf("taille %u (ncell=%d)\n", nelem, m->_nceli);
  
  double* Fx = readAndUncompress(ptrFile, nelem, elmin);
  double* Fy = readAndUncompress(ptrFile, nelem, elmin);
  double* Fz = readAndUncompress(ptrFile, nelem, elmin);
  
  char* nomFieldx = new char [STRINGLENGTH+1]; strcpy(nomFieldx, nomField); strcat(nomFieldx, "X");
  char* nomFieldy = new char [STRINGLENGTH+1]; strcpy(nomFieldy, nomField); strcat(nomFieldy, "Y");
  char* nomFieldz = new char [STRINGLENGTH+1]; strcpy(nomFieldz, nomField); strcat(nomFieldz, "Z");
  if ((int)nelem == m->_nsom) 
  { 
    m->_nfieldNames.push_back(nomFieldx); m->_nfields.push_back(Fx);
    m->_nfieldNames.push_back(nomFieldy); m->_nfields.push_back(Fy);
    m->_nfieldNames.push_back(nomFieldz); m->_nfields.push_back(Fz);  
  }
  else if ((int)nelem >= m->_nceli) 
  { 
    m->_cfieldNames.push_back(nomFieldx); m->_cfields.push_back(Fx);
    m->_cfieldNames.push_back(nomFieldy); m->_cfields.push_back(Fy);
    m->_cfieldNames.push_back(nomFieldz); m->_cfields.push_back(Fz);
  }
  else
  {
    delete [] Fx; delete [] Fy; delete [] Fz;
  }
}

// Lecture d'un tenseur volumique
void readTenseurVolumique(FILE* ptrFile, std::map<unsigned int, mesh*>& meshes)
{
  // champ moyen ou instantane
  unsigned char moyen;
  fread(&moyen, sizeof(unsigned char), 1, ptrFile); //moyen = UIBE(moyen);
  printf("Champ moyen=%hhu\n", moyen);
  
  unsigned int nd;
  fread(&nd, sizeof(unsigned int), 1, ptrFile); nd = UIBE(nd);
  double temps;
  fread(&temps, sizeof(double), 1, ptrFile); temps = DBE(temps);
  printf("temps=" SF_F_ "\n", temps);
  mesh* m = checkMesh(nd, meshes);

  unsigned char classe;
  fread(&classe, sizeof(unsigned char), 1, ptrFile); //classe = UIBE(classe);
  printf("Champ de classe=%hhu\n", classe);
  
  // Unknown thing that is here (not in spec)
  char toto[20];
  fread(toto, sizeof(char), 8, ptrFile);
  
  // nom du champ
  char nomField[STRINGLENGTH];
  E_Int i = 0; char c;
  while ((c = fgetc(ptrFile)) != '\0' && i < STRINGLENGTH-1) { nomField[i] = c; i++; }
  nomField[i] = '\0';
  
  unsigned int elmin; unsigned int nelem;
  readNElemElmin(ptrFile, nelem, elmin);
  
  printf("nom champ=%s\n", nomField);
  printf("taille %u (ncell=%u)\n", nelem, m->_nceli);
  double* Fxx = readAndUncompress(ptrFile, nelem, elmin);
  double* Fxy = readAndUncompress(ptrFile, nelem, elmin);
  double* Fxz = readAndUncompress(ptrFile, nelem, elmin);
  double* Fyx = readAndUncompress(ptrFile, nelem, elmin);
  double* Fyy = readAndUncompress(ptrFile, nelem, elmin);
  double* Fyz = readAndUncompress(ptrFile, nelem, elmin);
  double* Fzx = readAndUncompress(ptrFile, nelem, elmin);
  double* Fzy = readAndUncompress(ptrFile, nelem, elmin);
  double* Fzz = readAndUncompress(ptrFile, nelem, elmin);
  
  char* nomFieldxx = new char [STRINGLENGTH+2]; strcpy(nomFieldxx, nomField); strcat(nomFieldxx, "XX");
  char* nomFieldxy = new char [STRINGLENGTH+2]; strcpy(nomFieldxy, nomField); strcat(nomFieldxy, "XY");
  char* nomFieldxz = new char [STRINGLENGTH+2]; strcpy(nomFieldxz, nomField); strcat(nomFieldxz, "XZ");
  
  char* nomFieldyx = new char [STRINGLENGTH+2]; strcpy(nomFieldyx, nomField); strcat(nomFieldyx, "YX");
  char* nomFieldyy = new char [STRINGLENGTH+2]; strcpy(nomFieldyy, nomField); strcat(nomFieldyy, "YY");
  char* nomFieldyz = new char [STRINGLENGTH+2]; strcpy(nomFieldyz, nomField); strcat(nomFieldyz, "YZ");
  
  char* nomFieldzx = new char [STRINGLENGTH+2]; strcpy(nomFieldzx, nomField); strcat(nomFieldzx, "ZX");
  char* nomFieldzy = new char [STRINGLENGTH+2]; strcpy(nomFieldzy, nomField); strcat(nomFieldzy, "ZY");
  char* nomFieldzz = new char [STRINGLENGTH+2]; strcpy(nomFieldzz, nomField); strcat(nomFieldzz, "ZZ");
  
  if ((int)nelem == m->_nsom) 
  { 
    m->_nfieldNames.push_back(nomFieldxx); m->_nfields.push_back(Fxx);
    m->_nfieldNames.push_back(nomFieldxy); m->_nfields.push_back(Fxy);
    m->_nfieldNames.push_back(nomFieldxz); m->_nfields.push_back(Fxz);
    m->_nfieldNames.push_back(nomFieldyx); m->_nfields.push_back(Fyx);
    m->_nfieldNames.push_back(nomFieldyy); m->_nfields.push_back(Fyy);
    m->_nfieldNames.push_back(nomFieldyz); m->_nfields.push_back(Fyz);
    m->_nfieldNames.push_back(nomFieldzx); m->_nfields.push_back(Fzx);
    m->_nfieldNames.push_back(nomFieldzy); m->_nfields.push_back(Fzy);
    m->_nfieldNames.push_back(nomFieldzz); m->_nfields.push_back(Fzz);
  }
  else if ((int)nelem >= m->_nceli)
  { 
    m->_cfieldNames.push_back(nomFieldxx); m->_cfields.push_back(Fxx);
    m->_cfieldNames.push_back(nomFieldxy); m->_cfields.push_back(Fxy);
    m->_cfieldNames.push_back(nomFieldxz); m->_cfields.push_back(Fxz);
    m->_cfieldNames.push_back(nomFieldyx); m->_cfields.push_back(Fyx);
    m->_cfieldNames.push_back(nomFieldyy); m->_cfields.push_back(Fyy);
    m->_cfieldNames.push_back(nomFieldyz); m->_cfields.push_back(Fyz);
    m->_cfieldNames.push_back(nomFieldzx); m->_cfields.push_back(Fzx);
    m->_cfieldNames.push_back(nomFieldzy); m->_cfields.push_back(Fzy);
    m->_cfieldNames.push_back(nomFieldzz); m->_cfields.push_back(Fzz);
  }  
  else
  {
    delete [] Fxx; delete [] Fxy; delete [] Fxz;
    delete [] Fyx; delete [] Fyy; delete [] Fyz;
    delete [] Fzx; delete [] Fzy; delete [] Fzz;
  }
}

//=============================================================================
/* 
   arcread - read binary cedre format (archive)
   IN: file: file name,
   OUT: varString: variables string
   OUT: structField: field for each structured zones,
   OUT: ni, nj, nk: number of points of each structured zones,
   OUT: unstructField: field for each unstructured zones, 
   OUT: connectivity: connectivity for each unstructured zones,
   OUT: eltType: eltType for each unstructured zones.

   eltType is:
   1: BAR
   2: TRI
   3: QUAD
   4: TETRA
   7: HEXA
   8: NGON
   return 1 if failure.
   return 0 if ok.

*/
//=============================================================================
E_Int K_IO::GenIO::arcread(
  char* file, char*& varString,
  vector<FldArrayF*>& structField,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connectivity,
  vector<E_Int>& eltType, vector<char*>& zoneNames,
  char*& centerVarString,
  vector<FldArrayF*>& centerStructField,
  vector<FldArrayF*>& centerUnstructField
  )
{
  /* File Opening */
  FILE* ptrFile;
  ptrFile = fopen(file, "rb");

  if (ptrFile == NULL)
  {
    printf("Warning: arcread: cannot open file %s.\n", file);
    return 1;
  }

  //printf("Error: arcread: not implemented.\n");
  unsigned char version;
  char name[STRINGLENGTH];
  char fileName[STRINGLENGTH];
  char unitName[STRINGLENGTH];
  char nomUti[STRINGLENGTH];
  char titre[STRINGLENGTH];
  char date[STRINGLENGTH];
  char machine[STRINGLENGTH];
  double tempsZero;
  double epsilon;
  unsigned char solverType;
  unsigned char dimField;
  int numMel;
  char nomMel[STRINGLENGTH];
  int nespeces;
  char nomEspece[STRINGLENGTH*10];
  int numGrp;
  char nomGrp[STRINGLENGTH];
  int nelem;
  char nomScalar[STRINGLENGTH];
  unsigned char typeSca;
  char unit[STRINGLENGTH];
  unsigned int read;
  unsigned int numUti;
  unsigned int nGe, nGr, ne, nr;
  unsigned int* dthGe=NULL;
  unsigned int *dthGr=NULL;
  unsigned int *dthe=NULL;
  double *dthr=NULL;
  std::map<unsigned int,mesh*> meshes;
  
  readBlockName(ptrFile, name);
  readVersion(ptrFile, &version);
  
  while (! feof(ptrFile))
  {
    readBlockName(ptrFile, name);
    if (strcmp(name, "VERSION") == 0) readVersion(ptrFile, &version);
    else if (strcmp(name, "TITRE") == 0) readTitle(ptrFile, titre, date, machine, tempsZero, epsilon);
    else if (strcmp(name, "SANS TITRE") == 0) readTitle(ptrFile, NULL, date, machine, tempsZero, epsilon);
    else if (strcmp(name, "GLOBAL") == 0) readGlobal(ptrFile, solverType, dimField);
    else if (strcmp(name, "ESPECES") == 0) readEspece(ptrFile, numMel, nomMel, nespeces, nomEspece);
    else if (strcmp(name, "SCALAIRES") == 0) readScalar(ptrFile, numGrp, nomGrp, nelem, nomScalar, typeSca);
    else if (strcmp(name, "UNITE") == 0) readUnit(ptrFile, unitName, unit);
    else if (strcmp(name, "THERMODYNAMIQUE") == 0) readThermo(ptrFile, read, fileName, nGe, nGr, ne, nr, dthGe, dthGr, dthe, dthr);
    else if (strcmp(name, "DOMUTIL") == 0) readDomutil(ptrFile, numUti, nomUti);
    else if (strcmp(name, "STRUCTURE") == 0) readStructure(ptrFile, meshes);
    else if (strcmp(name, "NUMEROTATION") == 0) readNumerotation(ptrFile, meshes);
    else if (strcmp(name, "CONNEXION") == 0) readConnexion(ptrFile, meshes);
    else if (strcmp(name, "MAILLAGE") == 0) readMaillage(ptrFile, meshes);
    else if (strcmp(name, "CELLULES") == 0) readCellules(ptrFile, meshes);
    else if (strcmp(name, "GRAVITE CELLULE") == 0) readGraviteCellule(ptrFile, meshes);
    else if (strcmp(name, "GRAVITE FACE") == 0) readGraviteFace(ptrFile, meshes);
    else if (strcmp(name, "SURFACE") == 0) readSurface(ptrFile, meshes);
    else if (strcmp(name, "VALEUR VOLUMIQUE") == 0) readValeurVolumique(ptrFile, meshes);
    else if (strcmp(name, "VECTEUR VOLUMIQUE") == 0) readVecteurVolumique(ptrFile, meshes);
    else if (strcmp(name, "TENSEUR VOLUMIQUE") == 0) readTenseurVolumique(ptrFile, meshes);
    else if (strcmp(name, "FIN_DU_FICHIER") == 0) { printf("Fin du fichier\n"); break; }
    else { printf("Block pas encore implemente: %s\n", name); break; }
  }
  delete [] dthGe; delete [] dthGr; delete [] dthe; delete [] dthr;
  fclose(ptrFile);
  
  // dimensionne les varStrings supposees identiques pour tous les blocs
  E_Int lx = 0;
  for (std::map<unsigned int,mesh*>::iterator it = meshes.begin(); it != meshes.end(); it++)
  { 
    mesh& m = *meshes[it->first];
    E_Int l = 6;
    for (size_t p = 0; p < m._nfieldNames.size(); p++) l += strlen(m._nfieldNames[p])+2;
    lx = std::max(l, lx);
  }
  varString = new char [lx];
  lx = 0;
  for (std::map<unsigned int,mesh*>::iterator it = meshes.begin(); it != meshes.end(); it++)
  { 
    mesh& m = *meshes[it->first];
    E_Int l = 1;
    for (size_t p = 0; p < m._cfieldNames.size(); p++) l += strlen(m._cfieldNames[p])+2;
    lx = std::max(l, lx);
  }
  centerVarString = new char [lx];
  printf("lx=" SF_D_ "\n", lx);
  
  // transforme les structures de maillages en arrays
  for (std::map<unsigned int,mesh*>::iterator it = meshes.begin(); it != meshes.end(); it++)
  {
    mesh& m = *meshes[it->first];
    char* zoneName = new char [128];
   
    if (m._type == 0) // structure
    {
      //structField.push_back(f);
      //strcpy(zoneName, m._nom);
      //zoneNames.push_back(zoneName);
    }
    else if (m._type == 2) // non structure NGON
    {
      strcpy(zoneName, m._nom);
      zoneNames.push_back(zoneName);
      // recopie des coord + champs en noeuds
      E_Int np = m._nsom;
      E_Int nvars = 3 + m._nfields.size();
      FldArrayF* f = new FldArrayF(np, nvars);
      E_Float* fx = f->begin(1);
      E_Float* fy = f->begin(2);
      E_Float* fz = f->begin(3);
      for (E_Int p = 0; p < np; p++) fx[p] = m._x[p];
      for (E_Int p = 0; p < np; p++) fy[p] = m._y[p];
      for (E_Int p = 0; p < np; p++) fz[p] = m._z[p];
      //printf("champs aux noeuds=%zu\n", m._nfields.size());
      for (size_t j = 0; j < m._nfields.size(); j++)
      {
        E_Float* fp = f->begin(4+j);
        E_Float* pf = m._nfields[j];
        for (E_Int p = 0; p < np; p++) fp[p] = pf[p];
      }
      // update varString    
      strcpy(varString, "x,y,z");
      for (size_t p = 0; p < m._nfields.size(); p++) 
      { strcat(varString, ","); strcat(varString, m._nfieldNames[p]); }
      //printf("nodes: %s\n", varString);
      
      // recopie des champs au centre + vire les ghost cells + sentinelle 
      //printf("champs aux centres=%zu, nceli=" SF_D_ "\n", m._cfields.size(), m._nceli);
      FldArrayF* fc = NULL;
      nvars = m._cfields.size();
      if (nvars > 0) fc = new FldArrayF(m._nceli, nvars);
      for (size_t j = 0; j < m._cfields.size(); j++)
      {
        E_Float* fp = fc->begin(j+1);
        E_Float* pf = m._cfields[j];
        for (E_Int p = 0; p < m._nceli; p++) fp[p] = pf[p];
        //for (E_Int p = 0; p < m._nceli; p++) fp[p] = 0.;   
      }
    
      // update centerVarString
      if (nvars > 0) strcpy(centerVarString, m._cfieldNames[0]);
      for (size_t p = 1; p < m._cfields.size(); p++)
      { strcat(centerVarString, ","); strcat(centerVarString, m._cfieldNames[p]); }
      printf("centers: %s\n", centerVarString);
      
      // Calcul de NFACE
      E_Int nf = m._nfac;
      FldArrayI* cFE = new FldArrayI(nf, 2);
      E_Int* cFEp = cFE->begin();
      E_Int ind;
      // start 1 et 0=exterior pour nous + supprime les ghost cells
      for (E_Int p = 0; p < nf; p++)
      {
        ind = m._pe[2*p]+1;
        if (ind > m._nceli) ind = 0; // tag ghost as exterior
        cFEp[p] = ind;
        ind = m._pe[2*p+1]+1;
        if (ind > m._nceli) ind = 0; // tag ghost as exterior
        cFEp[p+nf] = ind;
      } 
      //printf("FE: ");
      //for (E_Int p = 0; p < 2*nf; p++) printf(" " SF_D_ " ", cFEp[p]);
      //printf("\n");
      FldArrayI cNFace; E_Int nelts;
      K_CONNECT::connectFE2NFace(*cFE, cNFace, nelts);
      printf("I found " SF_D_ " elements\n", nelts);
      //printf("NFACE: ");
      //for (E_Int p = 0; p < cNFace.getSize(); p++) printf(" " SF_D_ " ", cNFace[p]);
      //printf("\n");
      for (E_Int p = 0; p < cNFace.getSize(); p++) { if (cNFace[p] > nf || cNFace[p] < 1) printf("DANGER " SF_D_ "\n",cNFace[p]);}
      delete cFE;
      E_Int sizeNFACE = cNFace.getSize();
      
      // Calcul de NGON
      E_Int sizeNGON = 0;
      for (E_Int p = 0; p < nf; p++) sizeNGON += (m._nbpoints[p]+1);
      printf("sizeNGON=" SF_D_ "\n", sizeNGON);
      FldArrayI* cn = new FldArrayI(sizeNGON+sizeNFACE+4);
      E_Int* cnp = cn->begin();
      cnp[0] = nf;
      cnp[1] = sizeNGON;
      E_Int* ptr = &cnp[2];
      unsigned int* ptr2 = m._numpoints;
      E_Int npf;
      for (E_Int p = 0; p < nf; p++)
      {
        npf = m._nbpoints[p];
        ptr[0] = npf;
        for (E_Int j = 0; j < npf; j++) ptr[1+j] = ptr2[j]+1;  
        ptr += npf+1; ptr2 += npf;
      }
      //printf("NGON: ");
      //for (E_Int p = 0; p < sizeNGON; p++) printf(" " SF_D_ " ",cnp[2+p]);
      //printf("\n");
      for (E_Int p = 0; p < sizeNGON; p++) { if (cnp[2+p] > nf || cnp[2+p] < 1) printf("DANGER " SF_D_ "\n",cnp[2+p]);}
    
      ptr[0] = nelts;
      ptr[1] = sizeNFACE;
      ptr += 2;
      E_Int* cNFacep = cNFace.begin();
      for (E_Int p = 0; p < sizeNFACE; p++) ptr[p] = cNFacep[p];
        
      //printf("final:\n");
      //for (E_Int p = 0; p < sizeNFACE+sizeNGON+4; p++) printf(" " SF_D_ " ", cnp[p]);
      //printf("\n");
      unstructField.push_back(f); connectivity.push_back(cn); eltType.push_back(8);
      centerUnstructField.push_back(fc);
    }
  }
  
  // delete meshes
  for (std::map<unsigned int,mesh*>::iterator it = meshes.begin(); it != meshes.end(); it++) delete meshes[it->first];
  
  return 0;
}

//=============================================================================
/*
  This routine enables binary write of archive file.
*/
//=============================================================================
E_Int K_IO::GenIO::arcwrite(char* file, char* dataFmt, char* varString,
      std::vector<E_Int>& ni, std::vector<E_Int>& nj, std::vector<E_Int>& nk,
      std::vector <FldArrayF*>& structField,
      std::vector<FldArrayF*>& unstructField,
      std::vector<FldArrayI*>& connect,
      std::vector< vector<E_Int> >& eltTypes,
      std::vector<char*>& zoneNames)
{
  printf("Error: arcwrite: not implemented.\n");
  return 1;
}

