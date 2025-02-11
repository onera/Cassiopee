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
// Code for dealing with deactivatedZones
#include "CPlotState.h"
#include "String/kstring.h"

// Clear deactivatedZones list
void CPlotState::clearDeactivatedZones() 
{
  chain_int* c = deactivatedZones;
  if (c == NULL) return;
  E_Int size = 0; while (c->next != NULL) { c = c->next; size++; }
  chain_int** ptrs = new chain_int* [size];
  c = deactivatedZones;
  size = 0; while (c->next != NULL) { ptrs[size] = c; c = c->next; size++; }
  for (E_Int i = 0; i < size; i++) free(ptrs[i]);
  delete [] ptrs;
  deactivatedZones = NULL;
}

// insert a zone number in deactivatedZones list (unique)
void CPlotState::insertDeactivatedZones(E_Int i)
{
  chain_int* c = deactivatedZones;
  chain_int* cp = c;
  E_Int present = 0;
  while (c != NULL)
  {
    if (c->value == i) { present = 1; break; }
    cp = c;
    c = c->next;
  }
  if (present == 1) return;

  struct chain_int* ci;
  ci = (struct chain_int*)malloc(sizeof(struct chain_int));
  ci->value = i;
  ci->next = NULL;
  if (cp == NULL) deactivatedZones = ci;
  else cp->next = ci;
} 

// remove a zone number in deactivatedZones list
void CPlotState::removeDeactivatedZones(E_Int i)
{
  chain_int* c = deactivatedZones;
  chain_int* cp = NULL;
  chain_int* cn = c;
  while (c != NULL)
  {
    if (c->value == i)
    { 
      if (cp == NULL) deactivatedZones = c->next;
      else cp->next = c->next; 
      cn = c;
      c = c->next;
      free(cn);
    }
    else 
    {
      cp = c;
      c = c->next;
    }
  }
}

// print list
void CPlotState::printDeactivatedZones()
{
  printf("==========================\n");
  chain_int* c = deactivatedZones;
  while (c != NULL)
  {
    printf("deactivated " SF_D_ "\n", c->value);
    c = c->next;
  }
  fflush(stdout);
}
