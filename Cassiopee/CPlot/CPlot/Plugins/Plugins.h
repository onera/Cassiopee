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
#ifndef _CPLOT_PLUGINS_H_
#define _CPLOT_PLUGINS_H_

class Data;

#include "../StructZone.h"

/* Define plugin chain functions */
struct chain_function_void
{
    char functionName[MAXSTRINGLENGTH];
    void (*f)(Data* d);
    struct chain_function_void* next;
};
struct chain_function_void2
{
    char functionName[MAXSTRINGLENGTH];
    char extension[MAXSTRINGLENGTH];
    void (*f)(Data* d, char *, char*, E_Int, E_Int, E_Int);
    struct chain_function_void2* next;
};
struct chain_function_void3
{
    char functionName[MAXSTRINGLENGTH];
    void (*f)(Data* d, E_Int, E_Int, E_Int, E_Int);
    struct chain_function_void3* next;
};
struct chain_function_int
{
    char functionName[MAXSTRINGLENGTH];
    char varName[MAXSTRINGLENGTH];
    int (*f)(Data* d, E_Int, E_Int, E_Int);
    struct chain_function_int* next;
};
struct chain_function_double
{
    char functionName[MAXSTRINGLENGTH];
    char varName[MAXSTRINGLENGTH];
    void (*f)(Data* d, double, float*, float*, float*);
    struct chain_function_double* next;
};

/* Define plugins library */
typedef struct
{
    struct chain_function_double* colorMap;
    struct chain_function_double* zoneColorMap;
    struct chain_function_void* lookFor;
    struct chain_function_void* addAVariable;
    struct chain_function_int* blanking;
    struct chain_function_void* selectNextZone;
    struct chain_function_void* selectPreviousZone;
    struct chain_function_void3* mouseClick;
    struct chain_function_void3* mouseMultipleClick;
    struct chain_function_void3* mouseAccurateClick;
    struct chain_function_void3* mouseRightClick;
    struct chain_function_void2* screenDump;
} Plugins;

// Add a variable plugins
void addPressureVariable();
double compVolumeOfBBOfCell(E_Int* ind, double* coord, E_Int dim);
void compInterpCellVertices(E_Int ic, E_Int jc, E_Int kc,
                            E_Int im, E_Int jm, E_Int km,          
                            E_Int* indtab);
void compMaxDiagOfCell(E_Int dim, E_Int* ind, double* coord, double* diag);
void setVariableToZero(E_Int ni, E_Int nj, E_Int nk, double* f);

// Blanking plugins (blanking.cpp)
int blankCellN(Data* d, E_Int p1, E_Int blank, E_Int zone);
int blankCellNF(Data* d, E_Int p1, E_Int blank, E_Int zone);
int blankStatus(Data* d, E_Int p1, E_Int blank, E_Int zone);

// Colormaps plugins (colormap.c)
void colBlueToRed(Data* d, double f, float* r, float* g, float* b);
void colGreenToRed(Data* d, double f, float* r, float* g, float* b);
void col2RGB(Data* d, double f, float* r, float* g, float* b);
void col2HSV(Data* d, double f, float* r, float* g, float* b);
void diverging(Data* d, double f, float* r, float* g, float* b);
void col3RGB(Data* d, double f, float* r, float* g, float* b);
void col3HSV(Data* d, double f, float* r, float* g, float* b);
void colMRGB(Data* d, double f, float* r, float* g, float* b);
void colMHSV(Data* d, double f, float* r, float* g, float* b);

// Look-for plugins (lookfor.cpp)
void lookForActiveZone(Data* d);
void lookForMaxValue(Data* d);
void lookForMinValue(Data* d);

// Select plugins (select.cpp)
void selectNextZoneIncr(Data* d);
void selectPreviousZoneIncr(Data* d);

// Mouse plugins
void mouseClickSelect(Data* d, E_Int button, E_Int etat, E_Int x, E_Int y);
void mouseClickMultipleSelect(Data* d, E_Int button, E_Int etat, E_Int x, E_Int y);
void mouseClickAccurateSelect(Data* d, E_Int button, E_Int etat, E_Int x, E_Int y);
void mouseRightClickSelect(Data* d, E_Int button, E_Int etat, E_Int x, E_Int y);
void mouseClickTag(Data* d, E_Int button, E_Int etat, E_Int x, E_Int y);
void mouseClickPoint(Data* d, E_Int button, E_Int etat, E_Int x, E_Int y);

// Screen dump plugins
void writePPMFile(Data* d, char *filename, char* buffer, 
                  E_Int win_width, E_Int win_height, E_Int mode);
void writePNGFile(Data* d, char *filename, char* buffer, 
                  E_Int width, E_Int height, E_Int mode);
void writeDPNGFile(Data* d, char *filename, char* buffer, 
                   E_Int width, E_Int height, E_Int mode);
void writeMPEGFrame(Data* d, char *filename, char* buffer, 
                    E_Int width, E_Int height, E_Int mode);
#endif
