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
  
//#if PY_VERSION_HEX >= 0x03000000
//  PyObject* pName = PyUnicode_FromString("KCore.installPath");
//#else
//  PyObject* pName = PyString_FromString("KCore.installPath");
//#endif
  //pName = Py_BuildValue("s", "KCore");
  //PyObject* pModule = PyImport_Import(pName);
  //Py_DECREF(pName);
  PyObject* pModule = PyImport_ImportModule("KCore.installPath");
  if (pModule != NULL)
  {
    PyObject* dict = PyModule_GetDict(pModule);

#if PY_VERSION_HEX >= 0x03000000    
    PyObject* o = PyDict_GetItem(dict, PyUnicode_FromString("libPath"));
    char* retChar = PyBytes_AsString(PyUnicode_AsUTF8String(o));
#else
    PyObject* o = PyDict_GetItem(dict, PyString_FromString("libPath"));
    char* retChar = PyString_AsString(o);
#endif
    Py_DECREF(pModule);
    return retChar;
  }
  else
  {
    PyErr_Print();
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
  //PyObject* pModule = PyImport_ImportModule("KCore.installPath");
  char* answer = NULL;
  if (pModule != NULL)
  {
    PyObject* dict = PyModule_GetDict(pModule);
    PyObject* func = PyDict_GetItemString(dict, "expanduser");
    if (func != NULL)
    {
      if (PyCallable_Check(func))
      {
#if PY_VERSION_HEX >= 0x03000000
        PyObject* sarg = PyUnicode_FromString("~");
#else  
        PyObject* sarg = PyString_FromString("~");
#endif
        PyObject* args = PyTuple_New(1);
        PyTuple_SetItem(args, 0, sarg);
        PyObject* rslt = PyObject_CallObject(func, args);
        if (rslt)
        {
#if PY_VERSION_HEX >= 0x03000000
          answer = PyBytes_AsString(PyUnicode_AsUTF8String(rslt));
#else 
          answer = PyString_AsString(rslt);
#endif
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
case 2254672957ul:
      date = 24230;
      break;
case 3986841136ul:
      date = 24233;
      break;
case 1758617645ul:
      date = 24236;
      break;
case 2761675278ul:
      date = 24239;
      break;
case 2449009237ul:
      date = 24242;
      break;
case 3657815600ul:
      date = 24245;
      break;
case 3100257509ul:
      date = 24248;
      break;
case 2103064291ul:
      date = 24251;
      break;
case 1376352586ul:
      date = 24254;
      break;
case 1179146691ul:
      date = 24257;
      break;
case 2163608599ul:
      date = 24260;
      break;
case 3963382966ul:
      date = 24263;
      break;
case 2751157766ul:
      date = 24266;
      break;
case 1219021024ul:
      date = 24269;
      break;
case 3621769571ul:
      date = 24272;
      break;
case 3239084798ul:
      date = 24275;
      break;
case 950339345ul:
      date = 24278;
      break;
case 2754839453ul:
      date = 24281;
      break;
case 1779041510ul:
      date = 24284;
      break;
case 1713909127ul:
      date = 24287;
      break;
case 1199589661ul:
      date = 24290;
      break;
case 337036135ul:
      date = 24293;
      break;
case 3963932026ul:
      date = 24296;
      break;
case 3865105664ul:
      date = 24299;
      break;
case 4080651234ul:
      date = 24302;
      break;
case 1547831496ul:
      date = 24305;
      break;
case 1855614844ul:
      date = 24308;
      break;
case 1396686470ul:
      date = 24311;
      break;
case 899898501ul:
      date = 24314;
      break;
case 3463601153ul:
      date = 24317;
      break;
case 2823258481ul:
      date = 24320;
      break;
case 3725733040ul:
      date = 24323;
      break;
case 1255338440ul:
      date = 24326;
      break;
case 170480166ul:
      date = 24329;
      break;
case 2890088604ul:
      date = 24332;
      break;
case 1145767977ul:
      date = 24335;
      break;
case 1100990204ul:
      date = 24338;
      break;
case 3474696052ul:
      date = 24341;
      break;
case 1072188893ul:
      date = 24344;
      break;
case 2079764244ul:
      date = 24347;
      break;
case 845831706ul:
      date = 24350;
      break;
case 3178746817ul:
      date = 24353;
      break;
case 3151977157ul:
      date = 24356;
      break;
case 675767303ul:
      date = 24359;
      break;
case 707697266ul:
      date = 24362;
      break;
case 768695784ul:
      date = 24365;
      break;
case 2727540856ul:
      date = 24368;
      break;
case 2569506697ul:
      date = 24371;
      break;
case 2997661211ul:
      date = 24374;
      break;
case 3381480308ul:
      date = 24377;
      break;
case 1046703113ul:
      date = 24380;
      break;
case 1985381292ul:
      date = 24383;
      break;
case 1709091957ul:
      date = 24386;
      break;
case 3957986472ul:
      date = 24389;
      break;
case 3629715785ul:
      date = 24392;
      break;
case 1193731144ul:
      date = 24395;
      break;
case 2437183996ul:
      date = 24398;
      break;
case 1738739596ul:
      date = 24401;
      break;
case 970915666ul:
      date = 24404;
      break;
case 1088106531ul:
      date = 24407;
      break;
case 4017057064ul:
      date = 24410;
      break;
case 4158157754ul:
      date = 24413;
      break;
case 1791215400ul:
      date = 24416;
      break;
case 2035456850ul:
      date = 24419;
      break;
case 336733221ul:
      date = 24422;
      break;
case 388967350ul:
      date = 24425;
      break;
case 526987461ul:
      date = 24428;
      break;
case 385164736ul:
      date = 24431;
      break;
case 1649488268ul:
      date = 24434;
      break;
case 1090332171ul:
      date = 24437;
      break;
case 892324115ul:
      date = 24440;
      break;
case 3024636181ul:
      date = 24443;
      break;
case 1969270626ul:
      date = 24446;
      break;
case 3326107284ul:
      date = 24449;
      break;
case 96295306ul:
      date = 24452;
      break;
case 3540921322ul:
      date = 24455;
      break;
case 4151366737ul:
      date = 24458;
      break;
case 2735071844ul:
      date = 24461;
      break;
case 91557144ul:
      date = 24464;
      break;
case 3402003840ul:
      date = 24467;
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
case 4084544040ul:
      date = 24230;
      break;
case 2106516990ul:
      date = 24233;
      break;
case 1482970376ul:
      date = 24236;
      break;
case 1046917570ul:
      date = 24239;
      break;
case 427653939ul:
      date = 24242;
      break;
case 3580299127ul:
      date = 24245;
      break;
case 922751228ul:
      date = 24248;
      break;
case 506380060ul:
      date = 24251;
      break;
case 4256941563ul:
      date = 24254;
      break;
case 4059307800ul:
      date = 24257;
      break;
case 3802797574ul:
      date = 24260;
      break;
case 841476716ul:
      date = 24263;
      break;
case 1974181017ul:
      date = 24266;
      break;
case 1241426664ul:
      date = 24269;
      break;
case 2861030408ul:
      date = 24272;
      break;
case 743693413ul:
      date = 24275;
      break;
case 21361186ul:
      date = 24278;
      break;
case 4034083453ul:
      date = 24281;
      break;
case 1726666656ul:
      date = 24284;
      break;
case 130731596ul:
      date = 24287;
      break;
case 1608813071ul:
      date = 24290;
      break;
case 3386893805ul:
      date = 24293;
      break;
case 2355437127ul:
      date = 24296;
      break;
case 191389608ul:
      date = 24299;
      break;
case 3533569198ul:
      date = 24302;
      break;
case 3141701867ul:
      date = 24305;
      break;
case 818227294ul:
      date = 24308;
      break;
case 432119547ul:
      date = 24311;
      break;
case 3083320994ul:
      date = 24314;
      break;
case 1608625058ul:
      date = 24317;
      break;
case 1767657641ul:
      date = 24320;
      break;
case 1397872803ul:
      date = 24323;
      break;
case 692364113ul:
      date = 24326;
      break;
case 193317690ul:
      date = 24329;
      break;
case 1675336494ul:
      date = 24332;
      break;
case 2740721513ul:
      date = 24335;
      break;
case 1681299554ul:
      date = 24338;
      break;
case 3737685251ul:
      date = 24341;
      break;
case 3814699638ul:
      date = 24344;
      break;
case 2635212780ul:
      date = 24347;
      break;
case 926392485ul:
      date = 24350;
      break;
case 2679280283ul:
      date = 24353;
      break;
case 3638919909ul:
      date = 24356;
      break;
case 3193901232ul:
      date = 24359;
      break;
case 1741857774ul:
      date = 24362;
      break;
case 3407096520ul:
      date = 24365;
      break;
case 3822778520ul:
      date = 24368;
      break;
case 1912035075ul:
      date = 24371;
      break;
case 3974930892ul:
      date = 24374;
      break;
case 277408033ul:
      date = 24377;
      break;
case 2046410764ul:
      date = 24380;
      break;
case 2381398943ul:
      date = 24383;
      break;
case 969151852ul:
      date = 24386;
      break;
case 673919182ul:
      date = 24389;
      break;
case 2274461645ul:
      date = 24392;
      break;
case 171072813ul:
      date = 24395;
      break;
case 3599718604ul:
      date = 24398;
      break;
case 2184127812ul:
      date = 24401;
      break;
case 3144528003ul:
      date = 24404;
      break;
case 1782266907ul:
      date = 24407;
      break;
case 490853195ul:
      date = 24410;
      break;
case 1192370159ul:
      date = 24413;
      break;
case 537690118ul:
      date = 24416;
      break;
case 1975970810ul:
      date = 24419;
      break;
case 1567814066ul:
      date = 24422;
      break;
case 2349189138ul:
      date = 24425;
      break;
case 568385842ul:
      date = 24428;
      break;
case 205974113ul:
      date = 24431;
      break;
case 4120020828ul:
      date = 24434;
      break;
case 4155232644ul:
      date = 24437;
      break;
case 4014027385ul:
      date = 24440;
      break;
case 2354778489ul:
      date = 24443;
      break;
case 2230614425ul:
      date = 24446;
      break;
case 3792946052ul:
      date = 24449;
      break;
case 1184901854ul:
      date = 24452;
      break;
case 3799718041ul:
      date = 24455;
      break;
case 4233655746ul:
      date = 24458;
      break;
case 3873069052ul:
      date = 24461;
      break;
case 2533953021ul:
      date = 24464;
      break;
case 3410196593ul:
      date = 24467;
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
case 1400346063ul:
      date = 24230;
      break;
case 2327856378ul:
      date = 24233;
      break;
case 1206688637ul:
      date = 24236;
      break;
case 647538557ul:
      date = 24239;
      break;
case 549474514ul:
      date = 24242;
      break;
case 1498253805ul:
      date = 24245;
      break;
case 75256645ul:
      date = 24248;
      break;
case 1218665153ul:
      date = 24251;
      break;
case 220268284ul:
      date = 24254;
      break;
case 2658210734ul:
      date = 24257;
      break;
case 1347300181ul:
      date = 24260;
      break;
case 3271439248ul:
      date = 24263;
      break;
case 195666387ul:
      date = 24266;
      break;
case 1666371216ul:
      date = 24269;
      break;
case 1777578699ul:
      date = 24272;
      break;
case 3780344397ul:
      date = 24275;
      break;
case 754497142ul:
      date = 24278;
      break;
case 2834675356ul:
      date = 24281;
      break;
case 3208995009ul:
      date = 24284;
      break;
case 2188276730ul:
      date = 24287;
      break;
case 1237840676ul:
      date = 24290;
      break;
case 1616377477ul:
      date = 24293;
      break;
case 4130055595ul:
      date = 24296;
      break;
case 2445336028ul:
      date = 24299;
      break;
case 1488803982ul:
      date = 24302;
      break;
case 452108249ul:
      date = 24305;
      break;
case 1528822056ul:
      date = 24308;
      break;
case 3674561877ul:
      date = 24311;
      break;
case 1328412540ul:
      date = 24314;
      break;
case 3246607602ul:
      date = 24317;
      break;
case 4265300082ul:
      date = 24320;
      break;
case 2809202887ul:
      date = 24323;
      break;
case 4086657025ul:
      date = 24326;
      break;
case 2011764298ul:
      date = 24329;
      break;
case 796301624ul:
      date = 24332;
      break;
case 4261225884ul:
      date = 24335;
      break;
case 476570760ul:
      date = 24338;
      break;
case 1985460674ul:
      date = 24341;
      break;
case 1511296821ul:
      date = 24344;
      break;
case 3117203824ul:
      date = 24347;
      break;
case 3218014696ul:
      date = 24350;
      break;
case 1149573357ul:
      date = 24353;
      break;
case 532017181ul:
      date = 24356;
      break;
case 3624211811ul:
      date = 24359;
      break;
case 2526102742ul:
      date = 24362;
      break;
case 3770841973ul:
      date = 24365;
      break;
case 3077740353ul:
      date = 24368;
      break;
case 1956388564ul:
      date = 24371;
      break;
case 704448775ul:
      date = 24374;
      break;
case 3592293535ul:
      date = 24377;
      break;
case 738631969ul:
      date = 24380;
      break;
case 993449289ul:
      date = 24383;
      break;
case 121823892ul:
      date = 24386;
      break;
case 3331785628ul:
      date = 24389;
      break;
case 1420054439ul:
      date = 24392;
      break;
case 3707928449ul:
      date = 24395;
      break;
case 2390569833ul:
      date = 24398;
      break;
case 61091734ul:
      date = 24401;
      break;
case 3752503347ul:
      date = 24404;
      break;
case 3087714420ul:
      date = 24407;
      break;
case 3097056987ul:
      date = 24410;
      break;
case 1520538290ul:
      date = 24413;
      break;
case 834602760ul:
      date = 24416;
      break;
case 612115721ul:
      date = 24419;
      break;
case 1238821535ul:
      date = 24422;
      break;
case 2559360030ul:
      date = 24425;
      break;
case 3221433766ul:
      date = 24428;
      break;
case 2450638000ul:
      date = 24431;
      break;
case 85875150ul:
      date = 24434;
      break;
case 1327702636ul:
      date = 24437;
      break;
case 18656081ul:
      date = 24440;
      break;
case 1922237765ul:
      date = 24443;
      break;
case 2441412439ul:
      date = 24446;
      break;
case 449252119ul:
      date = 24449;
      break;
case 1494538414ul:
      date = 24452;
      break;
case 2198523518ul:
      date = 24455;
      break;
case 2587266295ul:
      date = 24458;
      break;
case 3173535471ul:
      date = 24461;
      break;
case 3204171512ul:
      date = 24464;
      break;
case 559932868ul:
      date = 24467;
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
int K_KCORE::activation(const char* name)
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
