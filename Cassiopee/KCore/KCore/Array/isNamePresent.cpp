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
# include <string.h>
# include "Array/Array.h"

//=============================================================================
// Is name present in string?
// IN: name: nom a rechercher dans myString
// IN: myString: la chaine recherchee.
// Retourne -1 si name n'existe pas dans myString.
// Return i: si name existe et est en position i. La position etant 
// determinee suivant les virgules.
//=============================================================================
E_Int K_ARRAY::isNamePresent(const char* name, const char* myString)
{
  E_Int c = 0;
  E_Int compt;
  E_Int lenName = strlen(name);
  E_Int lenString = strlen(myString);
  E_Int diffLen = lenString-lenName;
  char l1, l2;

  for (c = 0; c <= diffLen; c++)
  {
    compt = 0;
    for (E_Int i = 0; i < lenName; i++)
    { if (name[i] == myString[c+i]) compt++; }
    if (compt == lenName)
    {
      if (c == 0) l1 = ' ';
      else l1 = myString[c-1];
      if (c == diffLen) l2 = ' ';
      else l2 = myString[c+lenName];
      if ((l1 == ',' || l1 == ' ') && (l2 == ',' || l2 == ' '))
        break;
    }
  }

  if (c > diffLen) return -1;
  
  // Compte les virgules
  compt = 0;
  for (E_Int p = 0; p < c; p++)
  { if (myString[p] == ',') compt++; }
  return compt;
}

//=============================================================================
// is CoordinateX present in string?
// IN: string: la chaine recherchee.
// Retourne -1 si name n'existe pas dans string.
// Return i: si name existe et est en position i. La position etant 
// determinee suivant les virgules.
//=============================================================================
E_Int K_ARRAY::isCoordinateXPresent(const char* string)
{
  E_Int r = K_ARRAY::isNamePresent("CoordinateX", string);
  if (r == -1) r = K_ARRAY::isNamePresent("x", string);
  if (r == -1) r = K_ARRAY::isNamePresent("X", string);
  return r;
}
//=============================================================================
// is CoordinateY present in myString?
// IN: myString: la chaine recherchee.
// Retourne -1 si name n'existe pas dans string.
// Return i: si name existe et est en position i. La position etant 
// determinee suivant les virgules.
//=============================================================================
E_Int K_ARRAY::isCoordinateYPresent(const char* myString)
{
  E_Int r = K_ARRAY::isNamePresent("CoordinateY", myString);
  if (r == -1) r = K_ARRAY::isNamePresent("y", myString);
  if (r == -1) r = K_ARRAY::isNamePresent("Y", myString);
  return r;
}
//=============================================================================
// is CoordinateZ present in string?
// IN: myString: la chaine recherchee.
// Retourne -1 si name n'existe pas dans string.
// Return i: si name existe et est en position i. La position etant 
// determinee suivant les virgules.
//=============================================================================
E_Int K_ARRAY::isCoordinateZPresent(const char* myString)
{
  E_Int r = K_ARRAY::isNamePresent("CoordinateZ", myString);
  if (r == -1) r = K_ARRAY::isNamePresent("z", myString);
  if (r == -1) r = K_ARRAY::isNamePresent("Z", myString);
  return r;
}
//=============================================================================
// is CellNatureField1 (1:discretise, 2:interpole, 0:masque) present 
// in myString?
// IN: myString: la chaine recherchee.
// Retourne -1 si name n'existe pas dans myString.
// Return i: si name existe et est en position i. La position etant 
// determinee suivant les virgules.
//=============================================================================
E_Int K_ARRAY::isCellNatureField1Present(const char* myString)
{
  E_Int r = K_ARRAY::isNamePresent("cellN", myString);
  if (r == -1) r = K_ARRAY::isNamePresent("cellnf", myString);
  return r;
}
//=============================================================================
// is CellNatureField2 (1:discretise, 0:masque, other values possible) 
// present in myString?
// IN: myString: la chaine recherchee.
// Retourne -1 si name n'existe pas dans myString.
// Return i: si name existe et est en position i. La position etant 
// determinee suivant les virgules.
//=============================================================================
E_Int K_ARRAY::isCellNatureField2Present(const char* myString)
{
  E_Int r = K_ARRAY::isNamePresent("cellN", myString);
  if (r == -1) r = K_ARRAY::isNamePresent("cellnf", myString);
  if (r == -1) r = K_ARRAY::isNamePresent("cellNF", myString);
  if (r == -1) r = K_ARRAY::isNamePresent("status", myString);
  if (r == -1) r = K_ARRAY::isNamePresent("ichim", myString);
  return r;
}
//=============================================================================
// is Density present in myString?
// IN: myString: la chaine recherchee.
// Retourne -1 si name n'existe pas dans myString.
// Return i: si name existe et est en position i. La position etant 
// determinee suivant les virgules.
//=============================================================================
E_Int K_ARRAY::isDensityPresent(const char* myString)
{
  E_Int r = K_ARRAY::isNamePresent("Density", myString);
  if (r == -1) r = K_ARRAY::isNamePresent("ro", myString);
  //if (r == -1) r = K_ARRAY::isNamePresent("RHO", myString);
  return r;
}
//=============================================================================
// is MomentumX present in myString?
// IN: myString: la chaine recherchee.
// Retourne -1 si name n'existe pas dans myString.
// Return i: si name existe et est en position i. La position etant 
// determinee suivant les virgules.
//=============================================================================
E_Int K_ARRAY::isMomentumXPresent(const char* myString)
{
  E_Int r = K_ARRAY::isNamePresent("MomentumX", myString);
  if (r == -1) r = K_ARRAY::isNamePresent("rou", myString);
  if (r == -1) r = K_ARRAY::isNamePresent("rovx", myString);
  //if (r == -1) r = K_ARRAY::isNamePresent("RHOU", myString);
  return r;
}
//=============================================================================
// is MomentumY present in myString?
// IN: myString: la chaine recherchee.
// Retourne -1 si name n'existe pas dans myString.
// Return i: si name existe et est en position i. La position etant 
// determinee suivant les virgules.
//=============================================================================
E_Int K_ARRAY::isMomentumYPresent(const char* myString)
{
  E_Int r = K_ARRAY::isNamePresent("MomentumY", myString);
  if (r == -1) r = K_ARRAY::isNamePresent("rov", myString);
  if (r == -1) r = K_ARRAY::isNamePresent("rovy", myString);
  //if (r == -1) r = K_ARRAY::isNamePresent("RHOV", myString);
  return r;
}
//=============================================================================
// is MomentumZ present in myString?
// IN: myString: la chaine recherchee.
// Retourne -1 si name n'existe pas dans myString.
// Return i: si name existe et est en position i. La position etant 
// determinee suivant les virgules.
//=============================================================================
E_Int K_ARRAY::isMomentumZPresent(const char* myString)
{
  E_Int r = K_ARRAY::isNamePresent("MomentumZ", myString);
  if (r == -1) r = K_ARRAY::isNamePresent("row", myString);
  if (r == -1) r = K_ARRAY::isNamePresent("rovz", myString);
  //if (r == -1) r = K_ARRAY::isNamePresent("RHOW", myString);
  return r;
}

//=============================================================================
// is EnergyStagnationDensity present in myString?
// IN: myString: la chaine recherchee.
// Retourne -1 si name n'existe pas dans myString.
// Return i: si name existe et est en position i. La position etant 
// determinee suivant les virgules.
//=============================================================================
E_Int K_ARRAY::isEnergyStagnationDensityPresent(const char* myString)
{
  E_Int r = K_ARRAY::isNamePresent("EnergyStagnationDensity", myString);
  if (r == -1) r = K_ARRAY::isNamePresent("roE", myString);
  if (r == -1) r = K_ARRAY::isNamePresent("roe", myString);
  //if (r == -1) r = K_ARRAY::isNamePresent("RHOE", myString);
  return r;
}
//=============================================================================
// is VelocityX present in myString?
// IN: myString: la chaine recherchee.
// Retourne -1 si name n'existe pas dans myString.
// Return i: si name existe et est en position i. La position etant 
// determinee suivant les virgules.
//=============================================================================
E_Int K_ARRAY::isVelocityXPresent(const char* myString)
{
  E_Int r = K_ARRAY::isNamePresent("VelocityX", myString);
  if (r == -1) r = K_ARRAY::isNamePresent("vx", myString);
  if (r == -1) r = K_ARRAY::isNamePresent("u", myString);
  if (r == -1) r = K_ARRAY::isNamePresent("U", myString);
  return r;
}
//=============================================================================
// is VelocityY present in myString?
// IN: myString: la chaine recherchee.
// Retourne -1 si name n'existe pas dans myString.
// Return i: si name existe et est en position i. La position etant 
// determinee suivant les virgules.
//=============================================================================
E_Int K_ARRAY::isVelocityYPresent(const char* myString)
{
  E_Int r = K_ARRAY::isNamePresent("VelocityY", myString);
  if (r == -1) r = K_ARRAY::isNamePresent("vy", myString);
  if (r == -1) r = K_ARRAY::isNamePresent("v", myString);
  if (r == -1) r = K_ARRAY::isNamePresent("V", myString);
  return r;
}
//=============================================================================
// is VelocityZ present in myString?
// IN: myString: la chaine recherchee.
// Retourne -1 si name n'existe pas dans myString.
// Return i: si name existe et est en position i. La position etant 
// determinee suivant les virgules.
//=============================================================================
E_Int K_ARRAY::isVelocityZPresent(const char* myString)
{
  E_Int r = K_ARRAY::isNamePresent("VelocityZ", myString);
  if (r == -1) r = K_ARRAY::isNamePresent("vz", myString);
  if (r == -1) r = K_ARRAY::isNamePresent("w", myString);
  if (r == -1) r = K_ARRAY::isNamePresent("W", myString);
  return r;
}
//=============================================================================
// is pressure present in myString?
// IN: myString: la chaine recherchee.
// Retourne -1 si name n'existe pas dans myString.
// Return i: si name existe et est en position i. La position etant 
// determinee suivant les virgules.
//=============================================================================
E_Int K_ARRAY::isTemperaturePresent(const char* myString)
{
  E_Int r = K_ARRAY::isNamePresent("Temperature", myString);
  if (r == -1) r = K_ARRAY::isNamePresent("Temp", myString);
  if (r == -1) r = K_ARRAY::isNamePresent("T", myString);
  return r;
}
//=============================================================================
// is pressure present in myString?
// IN: myString: la chaine recherchee.
// Retourne -1 si name n'existe pas dans myString.
// Return i: si name existe et est en position i. La position etant 
// determinee suivant les virgules.
//=============================================================================
E_Int K_ARRAY::isPressurePresent(const char* myString)
{
  E_Int r = K_ARRAY::isNamePresent("Pressure", myString);
  if (r == -1) r = K_ARRAY::isNamePresent("p", myString);
  if (r == -1) r = K_ARRAY::isNamePresent("P", myString);
  return r;
}

//=============================================================================
// is time present in myString?
// IN: myString: la chaine recherchee.
// Retourne -1 si name n'existe pas dans myString.
// Return i: si name existe et est en position i. La position etant 
// determinee suivant les virgules.
//=============================================================================
E_Int K_ARRAY::isTimePresent(const char* myString)
{
  E_Int r = K_ARRAY::isNamePresent("Time", myString);
  if (r == -1) r = K_ARRAY::isNamePresent("t", myString);
  return r;
}
