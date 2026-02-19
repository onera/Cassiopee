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
// Selection code
if (ptrState->selectionStyle == 0)
{
  if (zonep->selected == 1 && zonep->active == 1)            
    glColor4f(color2[0], color2[1], color2[2], ptrState->alpha); 
  else if (zonep->selected == 0 && zonep->active == 1)  
    glColor4f(color1[0], color1[1], color1[2], ptrState->alpha); 
  else glColor4f(color1[0], color1[1], color1[2], 0.2); 
}
else if (ptrState->selectionStyle == 1)
{ 
  color2[0] = color1[0]; color2[1] = color1[1]; color2[2] = color1[2];
  glColor4f(color1[0], color1[1], color1[2], ptrState->alpha);     
  float v2;                                                       
  if (ptrState->selectedZone <= 0) v2 = ptrState->alpha;
  else 
  { 
    if (zonep->selected == 0) v2 = 0.12;
    else v2 = ptrState->alpha; 
  } 
  if (zonep->active == 1)
    glColor4f(color1[0], color1[1], color1[2], v2);
}
else if (ptrState->selectionStyle == 2)
{
  color2[0] = color1[0]; color2[1] = color1[1]; color2[2] = color1[2];
  glColor4f(color1[0], color1[1], color1[2], ptrState->alpha);
  float v1, v2, alpha;
  if (zonep->previouslySelected == 0) v1 = 0.12;
  else v1 = ptrState->alpha;
  if (ptrState->selectedZone <= 0) v2 = ptrState->alpha;
  else
  {
    if (zonep->selected == 0) v2 = 0.12;
    else v2 = ptrState->alpha;
  }
  alpha = (v1-v2)/5.*ptrState->ktimer+v2;
  if (zonep->active == 1)
    glColor4f(color1[0], color1[1], color1[2], alpha);
}
