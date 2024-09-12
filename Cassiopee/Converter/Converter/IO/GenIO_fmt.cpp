/*    
    Copyright 2013-2024 Onera.

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
// Misc formatted IO routines

# include "GenIO.h"
# include "Array/Array.h"
# include "String/kstring.h"
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include "Def/DefFunction.h"
#include <errno.h>
#include <limits.h>

using namespace std;
using namespace K_FLD;

//=============================================================================
/* Verifie si word match dans buf.
   Ignore les majuscules.
   Retourne 1 si oui, 0 sinon.
*/
//=============================================================================
E_Int K_IO::GenIO::matchInString(char* buf, const char* word)
{
  E_Int l = strlen(word);
  E_Int p = strlen(buf);
  E_Int i, j, c;
  for (j = p-l; j >= 0; j--)
  {
    c = 0;
    for (i = l-1; i >=0 ; i--)
    {
      if (toupper(buf[j+i]) == toupper(word[i])) c++;
    }
    if (c == l) return 1;
  }
  return 0;
}

//=============================================================================
/* Lit un mot.
   On appelle separateur le dernier caractere lu apres le mot.
   Retourne 3 si un mot a ete lu et ; est le separateur
   Retourne 2 si un mot a ete lu et blank ou tab est le separateur
   Retourne 1 si un mot a ete lu et \n ou \r est le separateur 
   Retourne 0 si la fin du fichier a ete atteinte mais qu'un mot a pu etre lu
   Retourne -1 si la fin du fichier a ete atteinte et qu'aucun mot n'a pu
   etre lu.
   Passe les commentaires contenant un #
   Note: buf doit etre dimensionne a BUFSIZE
*/
//=============================================================================
E_Int K_IO::GenIO::readWord(FILE* ptrFile, char* buf)
{
  char t;
  t = fgetc(ptrFile);
  E_Int c = 0;
  while (t != EOF && c < BUFSIZE && (c == 0 || t != ' ') 
         && (c == 0 || t != '\n') 
         && (c == 0 || t != '\r')
         && (c == 0 || t != ';'))
  {
    if (t == '#') // skip comment
    {
      while (t != EOF && t != '\n' && t != '\r') t = fgetc(ptrFile);
      if (t == EOF) { buf[c] = '\0' ; return 0; }
    }
    if (t != ' ' && t != '\n' && t != '\r' && t != ';') { buf[c] = t; c++; }
    t = fgetc(ptrFile);
  }
  buf[c] = '\0';
  //printf("%s \n", buf);
  if (t == EOF && c != 0) return 0;
  else if (t == EOF) return -1;
  else if (t == ' ') return 2;
  else if (t == '\t') return 2;
  else if (t == ';') return 3;
  else return 1;
}

//=============================================================================
/*
  Lit le ptrFile jusqu'a rencontrer keyword ou eof.
  keyword doit etre en majuscule.
  Ignore les majuscules.
  Retourne 1 si keyword a ete trouve.
  Retourne 0 si eof a ete atteint.
*/
//=============================================================================
E_Int K_IO::GenIO::readGivenKeyword(FILE* ptrFile, const char* keyword)
{
  E_Int i, a;
  E_Int l = strlen(keyword);
  char* word = new char [l+1];
  for (i = 0; i < l; i++) word[i] = ' ';
  word[l] = '\0';

  char t;
  t = fgetc(ptrFile);
  
  while (t != EOF)
  {
    for (a = 0; a < l-1; a++) word[a] = word[a+1];
    word[l-1] = toupper(t);
    if (strcmp(word, keyword) == 0)
    {
      //printf("%s\n", word);
      delete [] word;
      return 1;
    }
    t = fgetc(ptrFile);
  }

  delete [] word;

  if (t == EOF) return 0;
  else return 1;
}

//=============================================================================
/*
  Lit le ptrFile jusqu'a rencontrer keyword1 ou keyword2 ou eof.
  keyword1 et keyword2 doicent etre en majuscule.
  Ignore les majuscules.
  Retourne 1 si keyword1 a ete trouve en premier.
  Retourne 2 si keyword2 a ete trouve en premier.
  Retourne 0 si eof a ete atteint.
*/
//=============================================================================
E_Int K_IO::GenIO::readGivenKeyword(FILE* ptrFile, const char* keyword1, 
                                    const char* keyword2)
{
  E_Int i, a;
  E_Int l1 = strlen(keyword1);
  E_Int l2 = strlen(keyword2);
  
  char* word1 = new char [l1+1];
  char* word2 = new char [l2+1];
  
  for (i = 0; i < l1; i++) word1[i] = ' ';
  word1[l1] = '\0';
  for (i = 0; i < l2; i++) word2[i] = ' ';
  word2[l2] = '\0';

  char t;
  t = fgetc(ptrFile);
  
  while (t != EOF)
  {
    for (a = 0; a < l1-1; a++) word1[a] = word1[a+1];
    word1[l1-1] = toupper(t);
    for (a = 0; a < l2-1; a++) word2[a] = word2[a+1];
    word2[l2-1] = toupper(t);
  
    if (strcmp(word1, keyword1) == 0)
    {
      //printf("%s\n", word1);
      delete [] word1; delete [] word2;
      return 1;
    }

    if (strcmp(word2, keyword2) == 0)
    {
      //printf("%s\n", word2);
      delete [] word1; delete [] word2;
      return 2;
    }

    t = fgetc(ptrFile);
  }

  delete [] word1; delete [] word2;

  if (t == EOF) return 0;
  else return 1;
}

//=============================================================================
/*
  Lit le ptrFile jusqu'a rencontrer keyword1 ou keyword2 ou keyword3 ou eof.
  keyword1, keyword2, keyword3 doivent etre en majuscule.
  Ignore les majuscules.
  Retourne 1 si keyword1 a ete trouve en premier.
  Retourne 2 si keyword2 a ete trouve en premier.
  Retourne 3 si keyword3 a ete trouve en premier.
  Retourne 0 si eof a ete atteint.
*/
//=============================================================================
E_Int K_IO::GenIO::readGivenKeyword(FILE* ptrFile, const char* keyword1, 
                                    const char* keyword2, const char *keyword3)
{
  E_Int i, a;
  E_Int l1 = strlen(keyword1);
  E_Int l2 = strlen(keyword2);
  E_Int l3 = strlen(keyword3);
  
  char* word1 = new char [l1+1];
  char* word2 = new char [l2+1];
  char* word3 = new char [l3+1];
  
  for (i = 0; i < l1; i++) word1[i] = ' ';
  word1[l1] = '\0';
  for (i = 0; i < l2; i++) word2[i] = ' ';
  word2[l2] = '\0';
  for (i = 0; i < l3; i++) word3[i] = ' ';
  word2[l3] = '\0';

  char t;
  t = fgetc(ptrFile);
  
  while (t != EOF)
  {
    for (a = 0; a < l1-1; a++) word1[a] = word1[a+1];
    word1[l1-1] = toupper(t);
    for (a = 0; a < l2-1; a++) word2[a] = word2[a+1];
    word2[l2-1] = toupper(t);
    for (a = 0; a < l3-1; a++) word3[a] = word3[a+1];
    word3[l3-1] = toupper(t);
  
    if (strcmp(word1, keyword1) == 0)
    {
      //printf("%s\n", word1);
      delete [] word1; delete [] word2; delete [] word3;
      return 1;
    }

    if (strcmp(word2, keyword2) == 0)
    {
      //printf("%s\n", word2);
      delete [] word1; delete [] word2; delete [] word3;
      return 2;
    }

    if (strcmp(word3, keyword3) == 0)
    {
      //printf("%s\n", word3);
      delete [] word1; delete [] word2; delete [] word3;
      return 3;
    }

    t = fgetc(ptrFile);
  }

  delete [] word1; delete [] word2; delete [] word3;

  if (t == EOF) return 0;
  else return 1;
}

//=============================================================================
/*
  Lit le ptrFile jusqu'a recontrer '=' ou atteindre la taille max du buffer.
  En cas de succes, retourne 0 et le Keyword dans buf.
  Sinon retourne 1.
  Note: buf doit etre dimensionne a BUFSIZE
*/
//=============================================================================
E_Int K_IO::GenIO::readKeyword(FILE* ptrFile, char* buf)
{
  E_Int c = 0;
  char t;

  t = fgetc(ptrFile);

  while (c < BUFSIZE && t != '=' && t != EOF)
  {
    buf[c] = t;
    if (c>=4 && (t==' '||t==','||t=='\n') && strcmp(std::string(&buf[c-4], 4).c_str(),"ZONE")==0) break; // "ZONE" est un keyword non suivi de "="
    if (c>=4 && (t==' '||t==','||t=='\n') && strcmp(std::string(&buf[c-4], 4).c_str(),"TEXT")==0) break; // "TEXT" est un keyword non suivi de "="
    c++;
    t = fgetc(ptrFile);
  }
  
  if (c < BUFSIZE && t != EOF)
  {
    buf[c] = '\0';
    return 0;
  }
  else return 1;
}

//=============================================================================
/*
  Lit le ptrFile jusqu'a recontrer '=' ou atteindre la taille max du buffer.
  Essai de reconnaitre un keyword dans le buffer.
  Les knowKeywords doivent etre en majuscules.
  Ignore les majuscules.
  En cas de succes, retourne 0, prevData et keyword.
  En cas d'echec retourne 1.
*/
//=============================================================================
E_Int K_IO::GenIO::readDataAndKeyword(FILE* ptrFile, char* buf, 
                                      list<const char*>& knownKeywords,
                                      char* prevData, char* keyword)
{
  E_Int l, p, i, j;
  E_Int c = 0;

  readKwd:
  E_Int ret = readKeyword(ptrFile, buf);

  if (ret != 0)
  {
    buf[BUFSIZE] = '\0';
    strcpy(prevData, buf);
    return ret; // echec
  }

  // Find keyword in buf
  list<const char*>::iterator it;
  E_Int found = 0;
  for (it = knownKeywords.begin(); it != knownKeywords.end(); it++)
  {
    const char* k = *it;
    l = strlen(k);
    p = strlen(buf);
    for (j = p-l; j >= K_FUNC::E_max(0,p-l-2); j--)
    {
      c = 0;
      for (i = l-1; i >=0 ; i--)
      {
        if (toupper(buf[j+i]) == k[i]) c++;
      }
      if (c == l)
      {
        found = 1;
        strcpy(keyword, k); // stocke le keyword
        goto trouve;
      }
    }
  }

  if (found == 0) // pas un mot cle reconnu, continue a la recherche d'un nouveau mot cle
  {
    goto readKwd;
    //strcpy(prevData, buf);
    //return 1; // echec
  }

  // stocke la prevdata
  trouve:
  for (i = 0; i < j; i++) prevData[i] = buf[i];
  prevData[j] = '\0';

  return 0;
}

//=============================================================================
/*
  Lit un double dans le fichier.
  IN: ptrFile
  OUT: value
  IN: formatLength: si -1, pas utilise, sinon le format a formatLength 
  caracteres.
  Lit une valeur en caracteres, jusqu'a ce que ' ', '\n', ',', '<' ou '>' 
  ou tab est trouve ou formatLength est atteint.
  si eof est atteint, retourne 2
  si blank ou tab est le separateur, retourne 1
  si aucune valeur n'est lue retourne -1
  sinon retourne 0
  Some fortran have exponent as 'D', it is replaced with 'E'
*/
//=============================================================================
E_Int K_IO::GenIO::readDouble(FILE* ptrFile, E_Float& value, 
                              E_Int formatLength)
{
  E_Int i = 0;
  E_Int fc;
  char c;
  char number[256];

  if (formatLength == -1) // no formatLength imposed
  {
    c = 'a';
    while (c != EOF && (c != ' ' || i == 0) && (c != '\n' || i == 0) &&
           (c != '\r' || i == 0) && (c != ',' || i == 0) && 
           (c != '>' || i == 0) && (c != ';' || i == 0) &&
           (c != '<' || i == 0) && (c != '\v' || i == 0) && 
           (c != '\t' || i == 0) && (c != '(' || i == 0))
    {
      c = fgetc(ptrFile);
      if (c != ' ' && c != '\n' && c != '\r' && c != ',' && c != '>' 
          && c != '<' && c != '\v' && c != '\t' 
          && c != ';' && c != '(')
      {number[i] = c; if (c == 'D') number[i] = 'E'; i++;} 
    }
  }
  else
  {
    fc = 0; c = 'a';
    while (c != EOF && fc < formatLength)
    {
      c = fgetc(ptrFile);
      if (c != ' ' && c != '\n' && c != '\r')
      {number[i] = c; if (c == 'D') number[i] = 'E'; i++;}
      if (c != '\n' && c != '\r') fc++;
    }
  }

  number[i] = '\0';
  value = strtod(number, NULL);
  //printf("number %s " SF_F_ " " SF_D_ "\n", number, value, formatLength);

  if (i == 0) return -1;
  else if (c == EOF) return 2;
  else if (c == ' ') return 1;
  else if (c == '\t') return 1;
  else return 0;
}

//=============================================================================
/*
  Lit un double dans buf a partir de la position pos.
  Renvoie la nouvelle position.
  IN: buf: le buffer a lire
  IN: size: la taille du buffer
  IN/OUT: pos: la position de lecture
  OUT: value: le double lu
*/
//=============================================================================
E_Int K_IO::GenIO::readDouble(char* buf, E_Int size, E_Int& pos, 
                              E_Float& value)
{
  char number[256];
  int i = 0;
  int c;
  
  if (pos >= size) return -2;

  c = 'a';
  while (pos < size && (c != ' ' || i == 0) && (c != '\n' || i == 0) &&
         (c != ',' || i == 0) && (c != '>' || i == 0) && 
         (c != '/' || i == 0) && (c != ';' || i == 0) &&
         (c != '<' || i == 0) && (c != '\v' || i == 0) && 
         (c != '\t' || i == 0) && (c != '(' || i == 0))
    {
      c = buf[pos]; pos++;
      if (c != ' ' && c != '\n' && c != '\r' && c != ',' && c != '>' 
          && c != '<' && c != '\v' && c != '\t' && c != '/' 
          && c != ';' && c != '(')
      {number[i] = c; if (c == 'D') number[i] = 'E'; i++;}
    }

  number[i] = '\0';
  errno = 0;
  value = strtod(number, NULL);
  //printf("number: " SF_F_ " %s\n", value, number);

  if (pos < size) return 2;
  else if (c == ' ') return 1;
  else if (c == '\t') return 1;
  else if (errno == ERANGE) return -1;
  else return 0;
}

//=============================================================================
/*
  Read int in file.
  IN: ptrFile
  OUT: value
  IN: formatLength: if -1, not used, else the value has formatLength chars.
  Read a value from chars, until a ' ' or a '\n' or a ',' or a '<' 
  or a '>' or a tab or a / or a ';'
  is found or formatLength is reached.
  if end of file is reached, return 2
  if blank or tab is the separator, return 1
  if no valid value is found return -1
  if value is invalid, return -1
  else return 0
  Some fortran have exponent as 'D', it is replaced with 'E'
*/
//=============================================================================
E_Int K_IO::GenIO::readInt(FILE* ptrFile, E_Int& value, 
                           E_Int formatLength)
{
  E_Int i = 0; // nbre de caracteres lus
  E_Int fc;
  char c;
  char number[256];

  if (formatLength == -1) // no formatLength imposed
  {
    c = 'a';
    while (c != EOF && (c != ' ' || i == 0) && (c != '\n' || i == 0) &&
           (c != ',' || i == 0) && (c != '>' || i == 0) && 
           (c != '/' || i == 0) && (c != ';' || i == 0) &&
           (c != '<' || i == 0) && (c != '\v' || i == 0) && 
           (c != '\t' || i == 0) && (c != '(' || i == 0))
    {
      c = fgetc(ptrFile);
      if (c != ' ' && c != '\n' && c != '\r' && c != ',' && c != '>' 
          && c != '<' && c != '\v' && c != '\t' && c != '/' 
          && c != ';' && c!= '(')
      {number[i] = c; if (c == 'D') number[i] = 'E'; i++;}  
    }
  }
  else
  {
    fc = 0; c = 'a';
    while (c != EOF && fc < formatLength)
    {
      c = fgetc(ptrFile);
      if (c != ' ' && c != '\n'&& c != '\r')
      {number[i] = c; if (c == 'D') number[i] = 'E'; i++;}
      if (c != '\n' && c != '\r') fc++;
    }
  }

  number[i] = '\0';
  errno = 0;
  value = strtol(number, NULL, 0);
  //printf("number: " SF_D_ " string: %s errno=" SF_D_ "\n", value, number, errno);

  if (i == 0) return -1; // pas de caractere valide lu
  else if (c == EOF) return 2; // end of file
  else if (c == ' ') return 1; // separateur blank
  else if (c == '\t') return 1; // separateur tab
  else if (errno == ERANGE) return -1; // pas un entier
  else return 0;
}

//=============================================================================
/*
  Read int in file written in hexa.
  See readInt
*/
//=============================================================================
E_Int K_IO::GenIO::readHexaInt(FILE* ptrFile, E_Int& value, 
                               E_Int formatLength)
{
  E_Int i = 0; // nbre de caracteres lus
  E_Int fc;
  char c;
  char number[256];

  if (formatLength == -1) // no formatLength imposed
  {
    c = 'a';
    while (c != EOF && (c != ' ' || i == 0) && (c != '\n' || i == 0) &&
           (c != ',' || i == 0) && (c != '>' || i == 0) && 
           (c != '/' || i == 0) && (c != ';' || i == 0) &&
           (c != '<' || i == 0) && (c != '\v' || i == 0) && 
           (c != '\t' || i == 0) && (c != '(' || i == 0))
    {
      c = fgetc(ptrFile);
      if (c != ' ' && c != '\n' && c != '\r' && c != ',' && c != '>' 
          && c != '<' && c != '\v' && c != '\t' && c != '/' 
          && c != ';' && c != '(')
      {number[i] = c; if (c == 'D') number[i] = 'E'; i++;}  
    }
  }
  else
  {
    fc = 0; c = 'a';
    while (c != EOF && fc < formatLength)
    {
      c = fgetc(ptrFile);
      if (c != ' ' && c != '\n'&& c != '\r')
      {number[i] = c; if (c == 'D') number[i] = 'E'; i++;}
      if (c != '\n' && c != '\r') fc++;
    }
  }

  number[i] = '\0';
  errno = 0;
  value = strtol(number, NULL, 16);
  //printf("number: " SF_D_ " %s\n", value, number);

  if (i == 0) return -1; // pas de caractere valide lu
  else if (c == EOF) return 2; // end of file
  else if (c == ' ') return 1; // separateur blank
  else if (c == '\t') return 1; // separateur tab
  else if (errno == ERANGE) return -1; // pas un entier
  else return 0;
}

//=============================================================================
/*
  Lit un entier dans buf a partir de la position pos.
  Renvoie la nouvelle position.
  IN: buf: le buffer a lire
  IN: size: la taille du buffer
  IN/OUT: pos: la position de lecture
  OUT: value: l'entier lu
*/
//=============================================================================
E_Int K_IO::GenIO::readInt(char* buf, E_Int size, E_Int& pos, E_Int& value)
{
  char number[256];
  int i = 0; // nbre de caracteres valides lus
  int c;
  
  if (pos >= size) return -2;

  c = 'a';
  while (pos < size && (c != ' ' || i == 0) && (c != '\n' || i == 0) &&
         (c != ',' || i == 0) && (c != '>' || i == 0) && 
         (c != '/' || i == 0) && (c != ';' || i == 0) &&
         (c != '<' || i == 0) && (c != '\v' || i == 0) && 
         (c != '\t' || i == 0) && (c != '(' || i == 0))
    {
      c = buf[pos]; pos++;
      if (c != ' ' && c != '\n' && c != '\r' && c != ',' && c != '>' 
          && c != '<' && c != '\v' && c != '\t' && c != '/' 
          && c != ';' && c != '(')
      {number[i] = c; if (c == 'D') number[i] = 'E'; i++;}
    }

  number[i] = '\0';
  errno = 0;
  value = strtol(number, NULL, 0);
  //printf("number : " SF_D_ " %s " SF_D_ "\n", value, number, i);
  
  if (i == 0) return -1; // pas de caractere valide lu
  else if (errno == ERANGE) return -1; // pas entier
  else if (c == ' ') return 1; // sep est blank
  else if (c == '\t') return 1; // sep est tab
  else if (c == '\n') return 2; // sep en \n
  else if (pos < size) return 3; // pos est inferieure a size
  else return 0; // ok
}

//=============================================================================
/* 
   Lit un tuple 12/13/24 ou 12/13 ou 12.
   Retourne la premiere valeur.
   Retourne 2 si eof, 1 si arret sur un blanc et 0 sinon. 
 */
//=============================================================================
E_Int K_IO::GenIO::readIntTuple(FILE* ptrFile, E_Int& value)
{
  E_Int i = 0;
  char c;
  char number[256];
  
  c = fgetc(ptrFile);
  // Avance jusqu'a un tuple lisible
  while (c != EOF && (c == ' ' || c == '\n' || c == '\r'))
  {
    c = fgetc(ptrFile);
  }

  // Lit la premiere valeur du tuple
  while (c != EOF && c != ' ' && c != '/' && c != '\n' && c != '\r')
  {
    {number[i] = c; if (c == 'D') number[i] = 'E'; i++;}
    c = fgetc(ptrFile);
  }
  number[i] = '\0';
  
  // Lit les autres valeurs (si il y en a)
  while (c != EOF && c != ' ' && c != '\n' && c != '\r')
  {
    c = fgetc(ptrFile);
  }
  
  value = strtol(number, NULL, 0);

  if (c == EOF) return 2;
  else if (c == ' ') return 1;
  else if (c == '\t') return 1;
  else return 0;
}

//=============================================================================
/* 
   Lit un tuple 12/13/24 ou 12/13.
   Retourne les deux premieres valeurs.
   Retourne 2 si eof, 1 si arret sur un blanc et 0 sinon. 
 */
//=============================================================================
E_Int K_IO::GenIO::readIntTuple2(FILE* ptrFile, E_Int& value1, E_Int& value2)
{
  E_Int i = 0;
  char c;
  char number[256];
  
  c = fgetc(ptrFile);
  // Avance jusqu'a un tuple lisible
  while (c != EOF && (c == ' ' || c == '\n' || c == '\r'))
  {
    c = fgetc(ptrFile);
  }

  // Lit la premiere valeur du tuple
  while (c != EOF && c != ' ' && c != '/' && c != '\n' && c != '\r')
  {
    {number[i] = c; if (c == 'D') number[i] = 'E'; i++;}
    c = fgetc(ptrFile);
  }
  number[i] = '\0';

  value1 = strtol(number, NULL, 0);

  c = fgetc(ptrFile); // passe /
  i = 0;
  while (c != EOF && c != ' ' && c != '/' && c != '\n' && c != '\r')
  {
    {number[i] = c; if (c == 'D') number[i] = 'E'; i++;}
    c = fgetc(ptrFile);
  }
  number[i] = '\0';

  value2 = strtol(number, NULL, 0);

  // Lit les autres valeurs (si il y en a)
  while (c != EOF && c != ' ' && c != '\n' && c != '\r')
  {
    c = fgetc(ptrFile);
  }
  
  if (c == EOF) return 2;
  else if (c == ' ') return 1;
  else if (c == '\t') return 1;
  else return 0;
}

//=============================================================================
/* 
   Lit un tuple 12/13/24.
   Retourne les trois valeurs.
   Retourne 2 si eof, 1 si arret sur un blanc et 0 sinon. 
 */
//=============================================================================
E_Int K_IO::GenIO::readIntTuple3(FILE* ptrFile, E_Int& value1, E_Int& value2, E_Int& value3)
{
  E_Int i = 0;
  char c;
  char number[256];
  
  value1 = -1; value2 = -1; value3 = -1;

  c = fgetc(ptrFile);
  // Avance jusqu'a un tuple lisible
  while (c != EOF && (c == ' ' || c == '\n' || c == '\r'))
  {
    c = fgetc(ptrFile);
  }

  // Lit la premiere valeur du tuple
  while (c != EOF && c != ' ' && c != '/' && c != '\n' && c != '\r')
  {
    {number[i] = c; if (c == 'D') number[i] = 'E'; i++;}
    c = fgetc(ptrFile);
  }
  number[i] = '\0';

  value1 = strtol(number, NULL, 0);

  if (c == '/')
  { // essai pour lire la deuxieme valeur
    c = fgetc(ptrFile); // passe /
    i = 0;
    while (c != EOF && c != ' ' && c != '/' && c != '\n' && c != '\r')
    {
      {number[i] = c; if (c == 'D') number[i] = 'E'; i++;}
      c = fgetc(ptrFile);
    }
    number[i] = '\0';

    value2 = strtol(number, NULL, 0);

    if (c == '/')
    {
      c = fgetc(ptrFile); // passe /
      i = 0;
      while (c != EOF && c != ' ' && c != '/' && c != '\n' && c != '\r')
      {
        {number[i] = c; if (c == 'D') number[i] = 'E'; i++;}
        c = fgetc(ptrFile);
      }
      number[i] = '\0';

      value3 = strtol(number, NULL, 0);
    }
  }
  
  if (c == EOF) return 2;
  else if (c == ' ') return 1;
  else if (c == '\t') return 1;
  else return 0;
}

//=============================================================================
/* Skip comment.
   Skip a comment line in a file.
   Don't move ptrFile if the current line is not a comment.
   IN: commentChar: comment character
   IN: ptrFile: current ptr file
   OUT: ptrFile: modified.
   Retourne -1 si EOF,
   Retourn 1 si commentaire,
   Retourne 0 sinon.
*/
//=============================================================================
E_Int K_IO::GenIO::skipComment(FILE*& ptrFile, char commentChar)
{
  E_LONG pos = KFTELL(ptrFile);
  char c = fgetc(ptrFile);
  if (c == EOF)
  { KFSEEK(ptrFile, pos, SEEK_SET); return -1; }
  if (c == commentChar)
  {
    while (c != EOF && c != '\n' && c != '\r') c = fgetc(ptrFile);
    return 1;
  }
  KFSEEK(ptrFile, pos, SEEK_SET);
  return 0;
}

//=============================================================================
/* skip line
   Va jusqu'a EOF ou \n (passe)
   Retourne -1 si EOF
   sinon retourne 1.
*/
//=============================================================================
E_Int K_IO::GenIO::skipLine(FILE*& ptrFile)
{
  char c;
  while (1)
  {
    c = fgetc(ptrFile);
    if (c == EOF) return -1;
    if (c == '\n' || c == '\r') break;
  }
  return 1;
}

//=============================================================================
// Supprime les blancs (' ') et les \n dans la chaine buf
//=============================================================================
void K_IO::GenIO::compressString(char* buf)
{
  E_Int c = 0;
  E_Int i = 0;
  while (buf[c] != '\0')
  {
    if (buf[c] != ' ' && buf[c] != '\n' && buf[c] != '\r')
    {buf[i] = buf[c]; i++;}
    c++;
  }
  buf[i] = '\0';
}

//=============================================================================
// Lit jusqu'a rencontrer \n ou eof ou a depasser la taille fournie
// Retourne 0: ligne lue
// Retourne -1: eof
// Retourne 1: taille depassee
//=============================================================================
E_Int K_IO::GenIO::readline(FILE*& ptrFile, char* buf, E_Int size)
{
  E_Int i = 0;
  char c;
  while (i < size)
  {
    c = fgetc(ptrFile);
    if (c == EOF) return -1;
    buf[i] = c; i++;
    if (c == '\n') break;
  }
  if (i < size) {buf[i] = '\0'; return 0;}
  else return 1;
}

//=============================================================================
E_Int K_IO::GenIO::readTwoCoordinates(FILE* ptrFile, E_Float* pt)
{
  char* buf = new char[BUFSIZE+1];
  readWord(ptrFile, buf);
  E_Int ret = readTwoCoordinates(buf, ptrFile, pt);
  delete [] buf;
  return ret;
}

//=============================================================================
E_Int K_IO::GenIO::readTwoCoordinates(char* buf, FILE* ptrFile, E_Float* pt)
{
  vector<char*> numbers;
  K_ARRAY::extractVars(buf, numbers);
  E_Int numbersSize = numbers.size();
  E_Float value;
  for (E_Int i = 0; i < numbersSize; i++)
  { 
    value = strtod(numbers[i], NULL);
    delete [] numbers[i];
    pt[i] = value;
  }
  numbers.clear();
  
  if (numbersSize == 1)
  {
    readWord(ptrFile, buf);
    K_ARRAY::extractVars(buf, numbers);
    value = strtod(numbers[0], NULL);
    pt[1] = value;
    for (size_t i = 0; i < numbers.size(); i++) delete [] numbers[i];
  }
  
  E_Int lastChar = buf[strlen(buf)-1];
  if (lastChar == '\"') return 0;
  else return 1;
}

//=============================================================================
/* Write triangles into a formatted tecplot file. */
//=============================================================================
E_Int K_IO::GenIO::tpwriteTriangles(char* file, FldArrayF& field, 
                                    FldArrayI& connect, E_Boolean add)
{
  assert(connect.getNfld() == 3); // triangles
  E_Int nfld = field.getNfld();
  
  // Open file
  FILE* ptr_file = NULL;
  if (add == false)
  {
    ptr_file = fopen(file, "w");
    if (ptr_file == NULL) 
    {
      printf("Warning: tpwriteTriangles: I can't open file %s\n", file);
      exit(0);
    }
  
    // Title
    fprintf(ptr_file, "TITLE=\"triangles\"\n");
    switch (nfld)
    {
      case (3): // mesh
        fprintf(ptr_file, "VARIABLES=x y z\n");
        break;
        
      case (4): // mesh and cellN
        fprintf(ptr_file,"VARIABLES=x y z cellN\n");
        break;
        
      case (9): // euler
        fprintf(ptr_file, "VARIABLES=x y z ro rou rov row roE cellN\n");
        break;
        
      case (10):
        fprintf(ptr_file, "VARIABLES=x y z ro rou rov row roE rok cellN\n");
        break;
        
      case (11):
        fprintf(ptr_file, "VARIABLES=x y z ro rou rov row roE rok roeps cellN\n");
        break;
        
      default:
        printf("Number of variables unknown in tpwriteTriangles.\n");
        return 1;
    }
  }
  else
  {
    ptr_file = fopen(file,"a");
  }

  // Write just one zone
  E_Int npts = field.getSize();
  E_Int nbElts = connect.getSize();
  fprintf(ptr_file, "ZONE N=" SF_D_ ", E=" SF_D_ ", F=FEPOINT, ET=TRIANGLE\n",
          npts, nbElts);
  // Write points
  for (E_Int np = 0; np < npts; np++)
  {
    switch (nfld)
    {
      case (3): // mesh
        fprintf(ptr_file, "%2.6g %2.6g %2.6g\n",
                field(np,1), field(np,2), field(np,3));
        break;
        
      case (4): // mesh + cellN
        fprintf(ptr_file, "%2.6g %2.6g %2.6g %2.6g\n",
                field(np,1), field(np,2), field(np,3), field(np,4));
        break;
        
      case (9):
        fprintf(ptr_file,
                "%2.6g %2.6g %2.6g %2.6g %2.6g %2.6g %2.6g %2.6g %2.6g\n",
                field(np,1), field(np,2), field(np,3), field(np,4), field(np,5),
                field(np,6), field(np,7), field(np,8), field(np,9));
        break;
      case (10):
        fprintf(ptr_file, "%2.6g %2.6g %2.6g %2.6g %2.6g %2.6g %2.6g %2.6g %2.6g %2.6g\n",
                field(np,1), field(np,2), field(np,3), field(np,4), field(np,5),
                field(np,6), field(np,7), field(np,8), field(np,9), field(np,10));
        break;
      case (11):
        fprintf(ptr_file, "%2.6g %2.6g %2.6g %2.6g %2.6g %2.6g %2.6g %2.6g %2.6g %2.6g %2.6g\n",
                field(np,1), field(np,2), field(np,3), field(np,4), field(np,5),
                field(np,6), field(np,7), field(np,8), field(np,9), field(np,10), field(np,11));
        break;
        
      default:
        printf("Warning: tpwrite: Number of variables unknown in tpwriteTriangles.\n");
        return 1;
    }
  }
  
  // Write connectivity
  for (E_Int i = 0; i < nbElts; i++)
  { 
    fprintf(ptr_file, SF_D3_ "\n", connect(i,1), connect(i,2), connect(i,3));
  }

  fclose(ptr_file);
  return 0;
}

//=============================================================================
/* Write quads into a formatted tecplot file. */
//=============================================================================
E_Int K_IO::GenIO::tpwriteQuads(char* file, FldArrayF& field, 
                                FldArrayI& connect, E_Boolean add)
{
  assert(connect.getNfld() == 4); // quads
  E_Int nfld = field.getNfld();
  
  // Open file
  FILE* ptr_file = NULL;
  if (add == false)
  {
    ptr_file = fopen(file, "w");
    if (ptr_file == NULL) 
    {
      printf("Warning: tpwriteQuads: I can't open file %s.\n",file);
      exit(0);
    }
  
    // Title
    fprintf(ptr_file, "TITLE=\"quads\"\n");
    switch (nfld)
    {
      case (3): // mesh
        fprintf(ptr_file, "VARIABLES=x y z\n");
        break;
        
      case (4): // mesh and cellN
        fprintf(ptr_file, "VARIABLES=x y z cellN\n");
        break;
        
      case (9): // euler
        fprintf(ptr_file, "VARIABLES=x y z ro rou rov row roE cellN\n");
        break;
        
      case (10):
        fprintf(ptr_file, "VARIABLES=x y z ro rou rov row roE rok cellN\n");
        break;
        
      case (11):
        fprintf(ptr_file, "VARIABLES=x y z ro rou rov row roE rok roeps cellN\n");
        break;
        
      default:
        printf("Number of variables unknown in tpwriteQuads.\n");
        return 1;
    }
  }
  else
  {
    ptr_file = fopen(file, "a");
  }

  // Write just one zone
  E_Int npts = field.getSize();
  E_Int nbElts = connect.getSize();
  fprintf(ptr_file, "ZONE N=" SF_D_ ", E=" SF_D_ ", F=FEPOINT, ET=QUADRILATERAL\n",
          npts, nbElts);
  // Write points
  for (E_Int np = 0; np < npts; np++)
  {
    switch (nfld)
    {
      case (3): // mesh
        fprintf(ptr_file,"%2.6g %2.6g %2.6g\n",
                field(np,1), field(np,2), field(np,3));
        break;
        
      case (4): // mesh + cellN
        fprintf(ptr_file,"%2.6g %2.6g %2.6g %2.6g\n",
                field(np,1), field(np,2), field(np,3), field(np,4));
        break;
        
      case (9):
        fprintf(ptr_file,"%2.6g %2.6g %2.6g %2.6g %2.6g %2.6g %2.6g %2.6g %2.6g\n",
                field(np,1), field(np,2), field(np,3), field(np,4), field(np,5),
                field(np,6), field(np,7), field(np,8), field(np,9));
        break;
      case (10):
        fprintf(ptr_file,"%2.6g %2.6g %2.6g %2.6g %2.6g %2.6g %2.6g %2.6g %2.6g %2.6g\n",
                field(np,1), field(np,2), field(np,3), field(np,4), field(np,5),
                field(np,6), field(np,7), field(np,8), field(np,9), field(np,10));
        break;
      case (11):
        fprintf(ptr_file,"%2.6g %2.6g %2.6g %2.6g %2.6g %2.6g %2.6g %2.6g %2.6g %2.6g %2.6g\n",
                field(np,1), field(np,2), field(np,3), field(np,4), field(np,5),
                field(np,6), field(np,7), field(np,8), field(np,9), field(np,10), field(np,11));
        break;
        
      default:
        printf("Warning: tpwrite: number of variables unknown in tpwriteQuads.\n");
        return 1;
    }
  }
  
  // Write connectivity
  for (E_Int i = 0; i < nbElts; i++)
  { 
    fprintf(ptr_file, SF_D4_ "\n", connect(i,1),
            connect(i,2), connect(i,3), connect(i,4));
  }

  fclose(ptr_file);
  return 0;
}

//=============================================================================
E_Int K_IO::GenIO::tpwriteText(char* file, FldArrayF& field, FldArrayI& number)
{
  char numbers[30];
  // Open file	
  FILE* ptr_file = fopen(file, "w");
  if (ptr_file == NULL) 
  {
    printf("Warning: tpwriteText: I can't open file %s.\n",file);
    return 1;
  }

  E_Int nfld = field.getNfld();
  if (nfld != 3)
  {
    printf("Warning: tpwriteText: only x,y,z field is possible.");
    return 1;
  }
  
  // Title
  fprintf(ptr_file, "TITLE=\"text\"\n");
  
  E_Int npts = field.getSize();

  // write text number for each point
  for (E_Int np = 0; np < npts; np++)
  {
    sprintf(numbers, SF_D_, number[np]);
    fprintf(ptr_file, "TEXT\n");
    fprintf(ptr_file, "CS=GRID\n");
    fprintf(ptr_file, "C=BLACK\n");
    fprintf(ptr_file, "S=LOCAL\n"); 
    fprintf(ptr_file, "X=%2.6g,Y=%2.6g,Z=%2.6g\n",
            field(np,1), field(np,2), field(np,3));
    fprintf(ptr_file, "HU=GRID\n");
    fprintf(ptr_file, "LS=1 AN=LEFT\n");
    fprintf(ptr_file, "BXM=20 LT=0.1 BXO=BLACK BXF=WHITE\n");
    fprintf(ptr_file, "F=HELV-BOLD\n");
    fprintf(ptr_file, "H=8.47832201068E-05 A=0\n");
    fprintf(ptr_file, "MFC=\"\"");
    fprintf(ptr_file, "T=\"%s\"", numbers);
  }

  fclose(ptr_file);
  return 0;
}

//==============================================================================
// Converti une string en int.
// Retourne 0 si fail, sinon retourne l'entier
// Vire le debut s'il commence par 0 ou autre chose.
//==============================================================================
E_Int K_IO::GenIO::convertString2Int(char* str)
{
  // Get read of starting trash
  char* s = str;
  while (s[0] != '\0' && s[0] != '1' && s[0] != '+' && s[0] != '-' &&
         s[0] != '2' && s[0] != '3' && s[0] != '4' && s[0] != '5' &&
         s[0] != '6' && s[0] != '7' && s[0] != '8' && s[0] != '9') s++;

  return strtol(s, NULL, 0);
}
