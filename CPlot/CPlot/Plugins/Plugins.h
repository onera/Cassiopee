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
    void (*f)(Data* d, char *, char*,  int, int, int);
    struct chain_function_void2* next;
};
struct chain_function_void3
{
    char functionName[MAXSTRINGLENGTH];
    void (*f)(Data* d, int, int, int, int);
    struct chain_function_void3* next;
};
struct chain_function_int
{
    char functionName[MAXSTRINGLENGTH];
    char varName[MAXSTRINGLENGTH];
    int (*f)(Data* d, int, int, int);
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
double compVolumeOfBBOfCell( int* ind, double* coord, int dim);
void compInterpCellVertices(int ic, int jc, int kc,
                            int im, int jm, int km,          
                            int* indtab);
void compMaxDiagOfCell( int dim, int* ind, double* coord, double* diag );
void setVariableToZero(int ni, int nj, int nk, double* f);

// Blanking plugins (blanking.cpp)
int blankCellN(Data* d, int p1, int blank, int zone);
int blankCellNF(Data* d, int p1, int blank, int zone);
int blankStatus(Data* d, int p1, int blank, int zone);

// Colormaps plugins (colormap.c)
void colBlueToRed(Data* d, double f, float* r, float* g, float* b);
void colGreenToRed(Data* d, double f, float* r, float* g, float* b);
void colGrey(Data* d, double f, float* r, float* g, float* b);
void colGrey2(Data* d, double f, float* r, float* g, float* b);
void diverging(Data* d, double f, float* r, float* g, float* b);

// Look-for plugins (lookfor.cpp)
void lookForActiveZone(Data* d);
void lookForMaxValue(Data* d);
void lookForMinValue(Data* d);

// Select plugins (select.cpp)
void selectNextZoneIncr(Data* d);
void selectPreviousZoneIncr(Data* d);

// Mouse plugins
void mouseClickSelect(Data* d, int button, int etat, int x, int y);
void mouseClickMultipleSelect(Data* d, int button, int etat, int x, int y);
void mouseClickAccurateSelect(Data* d, int button, int etat, int x, int y);
void mouseRightClickSelect(Data* d, int button, int etat, int x, int y);
void mouseClickTag(Data* d, int button, int etat, int x, int y);
void mouseClickPoint(Data* d, int button, int etat, int x, int y);

// Screen dump plugins
void writePPMFile(Data* d, char *filename, char* buffer, 
                  int win_width, int win_height, int mode);
void writePNGFile(Data* d, char *filename, char* buffer, 
                  int width, int height, int mode);
void writeDPNGFile(Data* d, char *filename, char* buffer, 
                   int width, int height, int mode);
void writeMPEGFrame(Data* d, char *filename, char* buffer, 
                    int width, int height, int mode);
#endif
