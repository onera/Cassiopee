#include "kcore.h"
#include <time.h>

//=============================================================================
PyObject* K_KCORE::activation(PyObject* self, PyObject* args)
{
  char* name;
  PyArg_ParseTuple(args, "|s", &name);
  int date = activation(name);
  return Py_BuildValue("l", date);
}

//=============================================================================
// Get the install path from python
//=============================================================================
char* installPath()
{
  Py_Initialize();
  PyObject* pName = PyString_FromString("KCore.installPath");
  PyObject* pModule = PyImport_Import(pName);
  Py_DECREF(pName);
  if (pModule != NULL) 
  {
    PyObject* dict = PyModule_GetDict(pModule);
    PyObject* o = PyDict_GetItem(dict, PyString_FromString("libPath"));
    char* retChar = PyString_AsString(o);
    Py_DECREF(pModule);
    return retChar;
  }
  else
  {
    printf("Warning: module KCore.installPath can not be found.\n");
  }
  return NULL;
}

//=============================================================================
// Retourne le USER HOME - on utilise os.path.expanduser('~')
//=============================================================================
char* getHome()
{
  Py_Initialize();
  PyObject* pName = PyString_FromString("os.path");
  PyObject* pModule = PyImport_Import(pName);
  Py_DECREF(pName);
  char* answer = NULL;
  if (pModule != NULL)
  {
    PyObject* dict = PyModule_GetDict(pModule);
    PyObject* func = PyDict_GetItemString(dict, "expanduser");
    if (func != NULL)
    {
      if (PyCallable_Check(func))
      {
        PyObject* sarg = PyString_FromString("~");
        PyObject* args = PyTuple_New(1);
        PyTuple_SetItem(args, 0, sarg);
        PyObject* rslt = PyObject_CallObject(func, args);
        if (rslt)
        {
          answer = PyString_AsString(rslt);
          Py_DECREF(rslt);
        }
        Py_DECREF(sarg);
        Py_DECREF(args);
      }
    }
    Py_DECREF(pModule);
  }
  return answer;
}

//=============================================================================
E_LONG getDate()
{
  time_t rawtime;
  time(&rawtime);
  struct tm timeinfo = *localtime(&rawtime);

  int year = 1900 + timeinfo.tm_year;
  int month = timeinfo.tm_mon + 1;
  E_LONG date = year*12 + month;
  return date;
}
//=============================================================================
E_LONG htol(char* hash)
{
  unsigned int out[8];
  char hl[5];
  int j = 0;
  for (int k = 0; k < 8; k++)
  {
    for (int i = 0; i < 4; i++) { hl[i] = hash[j]; j++; }
    hl[4] = '\0';
    sscanf(hl, "%x", &out[k]);
  }
  E_LONG a = out[0] + out[1]*65535;
  E_LONG b = out[2] + out[3]*65535;
  E_LONG c = out[4] + out[5]*65535;
  E_LONG d = out[6] + out[7]*65535;
  return (a ^ b) ^ (c ^ d); 
}

//================================================================================
// Check a key
// Return 1: valid
// Return 0: invalid
//=================================================================================
int checkKeyA(char* key)
{
  E_LONG date;
  // Rip key si > 32 chars
  if (key[1] == ':')
  {
    for (E_Int i = 0; i < 32; i++) key[i] = key[i+2];
  }
  key[33] = '\0';
  //printf("key: %s\n", key);
  E_LONG myDate = getDate();
  E_LONG h = htol(key);
  //printf("tol: %ld\n", h);
  switch(h)
  {
case 2154702663ul:
      date = 24189;
      break;
case 1305957468ul:
      date = 24192;
      break;
case 3957104267ul:
      date = 24195;
      break;
case 1457299810ul:
      date = 24198;
      break;
case 509231139ul:
      date = 24201;
      break;
case 3108323998ul:
      date = 24204;
      break;
case 2364333103ul:
      date = 24207;
      break;
case 3247634524ul:
      date = 24210;
      break;
case 3072401952ul:
      date = 24213;
      break;
case 3459961138ul:
      date = 24216;
      break;
case 393870032ul:
      date = 24219;
      break;
case 1289993895ul:
      date = 24222;
      break;
case 1412226470ul:
      date = 24225;
      break;
case 1533863928ul:
      date = 24228;
      break;
case 1176951505ul:
      date = 24231;
      break;
case 574653943ul:
      date = 24234;
      break;
case 907005866ul:
      date = 24237;
      break;
case 2417627297ul:
      date = 24240;
      break;
case 2403978034ul:
      date = 24243;
      break;
case 3195827959ul:
      date = 24246;
      break;
case 1093221515ul:
      date = 24249;
      break;
case 4116558970ul:
      date = 24252;
      break;
case 2593845004ul:
      date = 24255;
      break;
case 2284155964ul:
      date = 24258;
      break;
case 626070987ul:
      date = 24261;
      break;
case 2762522645ul:
      date = 24264;
      break;
case 3341739770ul:
      date = 24267;
      break;
case 1040582695ul:
      date = 24270;
      break;
case 4206652715ul:
      date = 24273;
      break;
case 1166957892ul:
      date = 24276;
      break;
case 3098869956ul:
      date = 24279;
      break;
case 2744217969ul:
      date = 24282;
      break;
case 2813827895ul:
      date = 24285;
      break;
case 2777347949ul:
      date = 24288;
      break;
case 3509588958ul:
      date = 24291;
      break;
case 1124469994ul:
      date = 24294;
      break;
case 2225045217ul:
      date = 24297;
      break;
case 3324040815ul:
      date = 24300;
      break;
case 2686601226ul:
      date = 24303;
      break;
case 4173732699ul:
      date = 24306;
      break;
case 3535830428ul:
      date = 24309;
      break;
case 3210673783ul:
      date = 24312;
      break;
case 1812854888ul:
      date = 24315;
      break;
case 96621137ul:
      date = 24318;
      break;
case 4035859245ul:
      date = 24321;
      break;
case 1916332457ul:
      date = 24324;
      break;
case 3766369608ul:
      date = 24327;
      break;
case 3117315211ul:
      date = 24330;
      break;
case 309652015ul:
      date = 24333;
      break;
case 2130858270ul:
      date = 24336;
      break;
case 2381604303ul:
      date = 24339;
      break;
case 598413268ul:
      date = 24342;
      break;
case 702178807ul:
      date = 24345;
      break;
case 2712139184ul:
      date = 24348;
      break;
case 3817799065ul:
      date = 24351;
      break;
case 1791347490ul:
      date = 24354;
      break;
case 2128210142ul:
      date = 24357;
      break;
case 660041193ul:
      date = 24360;
      break;
case 3551635065ul:
      date = 24363;
      break;
case 2600594363ul:
      date = 24366;
      break;
case 1589631007ul:
      date = 24369;
      break;
case 2525365401ul:
      date = 24372;
      break;
case 2981993704ul:
      date = 24375;
      break;
case 2142448064ul:
      date = 24378;
      break;
case 1394453312ul:
      date = 24381;
      break;
case 631799887ul:
      date = 24384;
      break;
case 231714010ul:
      date = 24387;
      break;
case 2881646831ul:
      date = 24390;
      break;
case 557286536ul:
      date = 24393;
      break;
case 470585653ul:
      date = 24396;
      break;
case 1588146127ul:
      date = 24399;
      break;
case 2972966536ul:
      date = 24402;
      break;
case 797696072ul:
      date = 24405;
      break;
case 296456586ul:
      date = 24408;
      break;
case 2501733501ul:
      date = 24411;
      break;
case 508328102ul:
      date = 24414;
      break;
case 3281939631ul:
      date = 24417;
      break;
case 1399684697ul:
      date = 24420;
      break;
case 4200570855ul:
      date = 24423;
      break;
case 1820605996ul:
      date = 24426;
      break;

    default: 
      return 0;
  }
  if (myDate < date) return date;
  return 0;
}

int checkKey0(char* key)
{
  E_LONG date;
  // Rip key si > 32 chars
  if (key[1] == ':')
  {
    for (E_Int i = 0; i < 32; i++) key[i] = key[i+2];
  }
  key[33] = '\0';
  //printf("key: %s\n", key);
  E_LONG myDate = getDate();
  E_LONG h = htol(key);
  //printf("tol: %ld\n", h);
  switch(h)
  {
case 2154702663ul:
      date = 24189;
      break;
case 1305957468ul:
      date = 24192;
      break;
case 3957104267ul:
      date = 24195;
      break;
case 1457299810ul:
      date = 24198;
      break;
case 509231139ul:
      date = 24201;
      break;
case 3108323998ul:
      date = 24204;
      break;
case 2364333103ul:
      date = 24207;
      break;
case 3247634524ul:
      date = 24210;
      break;
case 3072401952ul:
      date = 24213;
      break;
case 3459961138ul:
      date = 24216;
      break;
case 393870032ul:
      date = 24219;
      break;
case 1289993895ul:
      date = 24222;
      break;
case 1412226470ul:
      date = 24225;
      break;
case 1533863928ul:
      date = 24228;
      break;
case 1176951505ul:
      date = 24231;
      break;
case 574653943ul:
      date = 24234;
      break;
case 907005866ul:
      date = 24237;
      break;
case 2417627297ul:
      date = 24240;
      break;
case 2403978034ul:
      date = 24243;
      break;
case 3195827959ul:
      date = 24246;
      break;
case 1093221515ul:
      date = 24249;
      break;
case 4116558970ul:
      date = 24252;
      break;
case 2593845004ul:
      date = 24255;
      break;
case 2284155964ul:
      date = 24258;
      break;
case 626070987ul:
      date = 24261;
      break;
case 2762522645ul:
      date = 24264;
      break;
case 3341739770ul:
      date = 24267;
      break;
case 1040582695ul:
      date = 24270;
      break;
case 4206652715ul:
      date = 24273;
      break;
case 1166957892ul:
      date = 24276;
      break;
case 3098869956ul:
      date = 24279;
      break;
case 2744217969ul:
      date = 24282;
      break;
case 2813827895ul:
      date = 24285;
      break;
case 2777347949ul:
      date = 24288;
      break;
case 3509588958ul:
      date = 24291;
      break;
case 1124469994ul:
      date = 24294;
      break;
case 2225045217ul:
      date = 24297;
      break;
case 3324040815ul:
      date = 24300;
      break;
case 2686601226ul:
      date = 24303;
      break;
case 4173732699ul:
      date = 24306;
      break;
case 3535830428ul:
      date = 24309;
      break;
case 3210673783ul:
      date = 24312;
      break;
case 1812854888ul:
      date = 24315;
      break;
case 96621137ul:
      date = 24318;
      break;
case 4035859245ul:
      date = 24321;
      break;
case 1916332457ul:
      date = 24324;
      break;
case 3766369608ul:
      date = 24327;
      break;
case 3117315211ul:
      date = 24330;
      break;
case 309652015ul:
      date = 24333;
      break;
case 2130858270ul:
      date = 24336;
      break;
case 2381604303ul:
      date = 24339;
      break;
case 598413268ul:
      date = 24342;
      break;
case 702178807ul:
      date = 24345;
      break;
case 2712139184ul:
      date = 24348;
      break;
case 3817799065ul:
      date = 24351;
      break;
case 1791347490ul:
      date = 24354;
      break;
case 2128210142ul:
      date = 24357;
      break;
case 660041193ul:
      date = 24360;
      break;
case 3551635065ul:
      date = 24363;
      break;
case 2600594363ul:
      date = 24366;
      break;
case 1589631007ul:
      date = 24369;
      break;
case 2525365401ul:
      date = 24372;
      break;
case 2981993704ul:
      date = 24375;
      break;
case 2142448064ul:
      date = 24378;
      break;
case 1394453312ul:
      date = 24381;
      break;
case 631799887ul:
      date = 24384;
      break;
case 231714010ul:
      date = 24387;
      break;
case 2881646831ul:
      date = 24390;
      break;
case 557286536ul:
      date = 24393;
      break;
case 470585653ul:
      date = 24396;
      break;
case 1588146127ul:
      date = 24399;
      break;
case 2972966536ul:
      date = 24402;
      break;
case 797696072ul:
      date = 24405;
      break;
case 296456586ul:
      date = 24408;
      break;
case 2501733501ul:
      date = 24411;
      break;
case 508328102ul:
      date = 24414;
      break;
case 3281939631ul:
      date = 24417;
      break;
case 1399684697ul:
      date = 24420;
      break;
case 4200570855ul:
      date = 24423;
      break;
case 1820605996ul:
      date = 24426;
      break;

    default: 
      return 0;
  }
  if (myDate < date) return date;
  return 0;
}

int checkKey1(char* key)
{
  E_LONG date;
  // Rip key si > 32 chars
  if (key[1] == ':')
  {
    for (E_Int i = 0; i < 32; i++) key[i] = key[i+2];
  }
  key[33] = '\0';
  //printf("key: %s\n", key);
  E_LONG myDate = getDate();
  E_LONG h = htol(key);
  //printf("tol: %ld\n", h);
  switch(h)
  {
case 2154702663ul:
      date = 24189;
      break;
case 1305957468ul:
      date = 24192;
      break;
case 3957104267ul:
      date = 24195;
      break;
case 1457299810ul:
      date = 24198;
      break;
case 509231139ul:
      date = 24201;
      break;
case 3108323998ul:
      date = 24204;
      break;
case 2364333103ul:
      date = 24207;
      break;
case 3247634524ul:
      date = 24210;
      break;
case 3072401952ul:
      date = 24213;
      break;
case 3459961138ul:
      date = 24216;
      break;
case 393870032ul:
      date = 24219;
      break;
case 1289993895ul:
      date = 24222;
      break;
case 1412226470ul:
      date = 24225;
      break;
case 1533863928ul:
      date = 24228;
      break;
case 1176951505ul:
      date = 24231;
      break;
case 574653943ul:
      date = 24234;
      break;
case 907005866ul:
      date = 24237;
      break;
case 2417627297ul:
      date = 24240;
      break;
case 2403978034ul:
      date = 24243;
      break;
case 3195827959ul:
      date = 24246;
      break;
case 1093221515ul:
      date = 24249;
      break;
case 4116558970ul:
      date = 24252;
      break;
case 2593845004ul:
      date = 24255;
      break;
case 2284155964ul:
      date = 24258;
      break;
case 626070987ul:
      date = 24261;
      break;
case 2762522645ul:
      date = 24264;
      break;
case 3341739770ul:
      date = 24267;
      break;
case 1040582695ul:
      date = 24270;
      break;
case 4206652715ul:
      date = 24273;
      break;
case 1166957892ul:
      date = 24276;
      break;
case 3098869956ul:
      date = 24279;
      break;
case 2744217969ul:
      date = 24282;
      break;
case 2813827895ul:
      date = 24285;
      break;
case 2777347949ul:
      date = 24288;
      break;
case 3509588958ul:
      date = 24291;
      break;
case 1124469994ul:
      date = 24294;
      break;
case 2225045217ul:
      date = 24297;
      break;
case 3324040815ul:
      date = 24300;
      break;
case 2686601226ul:
      date = 24303;
      break;
case 4173732699ul:
      date = 24306;
      break;
case 3535830428ul:
      date = 24309;
      break;
case 3210673783ul:
      date = 24312;
      break;
case 1812854888ul:
      date = 24315;
      break;
case 96621137ul:
      date = 24318;
      break;
case 4035859245ul:
      date = 24321;
      break;
case 1916332457ul:
      date = 24324;
      break;
case 3766369608ul:
      date = 24327;
      break;
case 3117315211ul:
      date = 24330;
      break;
case 309652015ul:
      date = 24333;
      break;
case 2130858270ul:
      date = 24336;
      break;
case 2381604303ul:
      date = 24339;
      break;
case 598413268ul:
      date = 24342;
      break;
case 702178807ul:
      date = 24345;
      break;
case 2712139184ul:
      date = 24348;
      break;
case 3817799065ul:
      date = 24351;
      break;
case 1791347490ul:
      date = 24354;
      break;
case 2128210142ul:
      date = 24357;
      break;
case 660041193ul:
      date = 24360;
      break;
case 3551635065ul:
      date = 24363;
      break;
case 2600594363ul:
      date = 24366;
      break;
case 1589631007ul:
      date = 24369;
      break;
case 2525365401ul:
      date = 24372;
      break;
case 2981993704ul:
      date = 24375;
      break;
case 2142448064ul:
      date = 24378;
      break;
case 1394453312ul:
      date = 24381;
      break;
case 631799887ul:
      date = 24384;
      break;
case 231714010ul:
      date = 24387;
      break;
case 2881646831ul:
      date = 24390;
      break;
case 557286536ul:
      date = 24393;
      break;
case 470585653ul:
      date = 24396;
      break;
case 1588146127ul:
      date = 24399;
      break;
case 2972966536ul:
      date = 24402;
      break;
case 797696072ul:
      date = 24405;
      break;
case 296456586ul:
      date = 24408;
      break;
case 2501733501ul:
      date = 24411;
      break;
case 508328102ul:
      date = 24414;
      break;
case 3281939631ul:
      date = 24417;
      break;
case 1399684697ul:
      date = 24420;
      break;
case 4200570855ul:
      date = 24423;
      break;
case 1820605996ul:
      date = 24426;
      break;

    default: 
      return 0;
  }
  if (myDate < date) return date;
  return 0;
}

//=============================================================================
// Retourne 1 si la cle d'activation correspondant a name est valide
// Retourne 0 sinon
// La cle doit etre dans .CassiopeeKey
//=============================================================================
int K_KCORE::activation(char* name)
{
  char filePath[256];
  int ret;
  char key[128];
  FILE* f=NULL;

  // try top open home/.CassiopeeKey
  char* home = getHome();
  if (home != NULL)
  {
    strcpy(filePath, home);
    strcat(filePath, "/.CassiopeeKey");
    f = fopen(filePath, "r");
  }
  if (f != NULL)
  {
    while (!feof(f))
    {
      fgets(key, 127, f);
      //printf("read key home: %s\n", key);
      if (name == NULL && strlen(key) == 32) { ret = checkKey0(key); if (ret != 0) {fclose(f); return ret;}}
      else if (name == NULL && strncmp(key,"A",1) == 0) { ret = checkKeyA(key); if (ret != 0) {fclose(f); return ret;}}
      else if (name != NULL && strncmp(key,"A",1) == 0) { ret = checkKeyA(key); if (ret != 0) {fclose(f); return ret;}}
      else if (name != NULL && strncmp(name,"0",1) == 0 && strncmp(key,"0",1) == 0) { ret = checkKey0(key); if (ret != 0) {fclose(f); return ret;}}
      else if (name != NULL && strncmp(name,"1",1) == 0 && strncmp(key,"1",1) == 0) { ret = checkKey1(key); if (ret != 0) {fclose(f); return ret;}}
    }
    fclose(f);
  }

  // try to open installPath/.CassiopeeKey
  char* path = installPath();
  if (path == NULL) return 0; // FAIL
  strcpy(filePath, path);
  strcat(filePath, "/.CassiopeeKey");
  f = fopen(filePath, "r");
  
  if (f == NULL) return 0; // no key file found
  
  while (!feof(f))
  {
    fgets(key, 127, f);
    //printf("read key dist: %s, name=%s\n", key,name);
    if (name == NULL && strlen(key) == 32) { ret = checkKey0(key); if (ret != 0) {fclose(f); return ret;}}
    else if (name == NULL && strncmp(key,"A",1) == 0) { ret = checkKeyA(key); if (ret != 0) {fclose(f); return ret;}}
    else if (name != NULL && strncmp(key,"A",1) == 0) { ret = checkKeyA(key); if (ret != 0) {fclose(f); return ret;}}
    else if (name != NULL && strncmp(name,"0",1) == 0 && strncmp(key,"0",1) == 0) { ret = checkKey0(key); if (ret != 0) {fclose(f); return ret;}}
    else if (name != NULL && strncmp(name,"1",1) == 0 && strncmp(key,"1",1) == 0) { ret = checkKey1(key); if (ret != 0) {fclose(f); return ret;}}
  }
  fclose(f);
  return 0; 
}
