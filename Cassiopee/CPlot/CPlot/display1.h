  // Common part to display methods  
  PyObject* arrays;
  int dim;
  PyObject* modeObject;
  PyObject* scalarFieldObject;
  PyObject* vectorFieldObject1, *vectorFieldObject2, *vectorFieldObject3;
  int winx, winy;
  int displayBB, displayInfo, displayIsoLegend;
  int meshStyle, solidStyle, scalarStyle, vectorStyle, colormap, niso;
  char* colormapC1; char* colormapC2; char* colormapC3; PyObject* colormapC;
  E_Float xcam, ycam, zcam, xeye, yeye, zeye, dirx, diry, dirz, isoEdges;
  E_Float stereoDist, viewAngle, vectorScale, vectorDensity;
  int vectorNormalize, vectorShowSurface, vectorShape, vectorProjection;
  char* exportFile; char* exportResolution; int exportAA;
  PyObject* zoneNamesObject;
  PyObject* renderTagsObject;
  PyObject* isoScales;
  char* backgroundFile;
  int bgColor, shadow, dof, offscreen, stereo, frameBuffer, panorama;
  E_Float lightOffsetX, lightOffsetY;
  E_Float dofPower, gamma; int toneMapping;
  PyObject* posCamList; PyObject* posEyeList; PyObject* dirCamList; // only for ODS
  if (!PyArg_ParseTuple(args, "OiOOOOOiiiiiiiddiiiiisssOidO(ii)(ddd)(ddd)(ddd)disi(dd)iddiidissiOOiiOOO",
                        &arrays, &dim, &modeObject, &scalarFieldObject,
                        &vectorFieldObject1, &vectorFieldObject2, &vectorFieldObject3,
                        &displayBB, &displayInfo, &displayIsoLegend,
                        &meshStyle, &solidStyle, &scalarStyle,
                        &vectorStyle, &vectorScale, &vectorDensity, &vectorNormalize, 
                        &vectorShowSurface, &vectorShape, &vectorProjection, 
                        &colormap, &colormapC1, &colormapC2, &colormapC3, &colormapC,
                        &niso, &isoEdges, &isoScales,
                        &winx, &winy, 
                        &xcam, &ycam, &zcam,
                        &xeye, &yeye, &zeye,
                        &dirx, &diry, &dirz, &viewAngle, 
                        &bgColor, &backgroundFile,
                        &shadow, &lightOffsetX, &lightOffsetY,
                        &dof, &dofPower, &gamma, &toneMapping, 
                        &stereo, &stereoDist, &panorama,
                        &exportFile, &exportResolution, &exportAA,
                        &zoneNamesObject, &renderTagsObject,
                        &frameBuffer, &offscreen,
                        &posCamList, &posEyeList, &dirCamList))
  {
    return NULL;
  }

  // Recuperation des noms de zones (eventuellement)
  vector<char*> zoneNames;
  getStringsFromPyObj(zoneNamesObject, zoneNames);

  // Recuperation des tags de render (eventuellement)
  vector<char*> renderTags;
  getStringsFromPyObj(renderTagsObject, renderTags);

  // Lecture des arrays
  vector<E_Int> res;
  vector<char*> structVarString; vector<char*> unstrVarString;
  vector<FldArrayF*> structF; vector<FldArrayF*> unstrF;
  vector<E_Int> nit; vector<E_Int> njt; vector<E_Int> nkt;
  vector<FldArrayI*> cnt;
  vector<char*> eltType;
  vector<PyObject*> objs, obju;

  char* varStringl; char* eltTypel;
  E_Int nil, njl, nkl, resl;
  FldArrayF* fl; FldArrayI* cnl;
  PyObject* o;
  E_Int nz = PyList_Size(arrays);
  for (E_Int i = 0; i < nz; i++)
  {
    o = PyList_GetItem(arrays, i);
    resl = K_ARRAY::getFromArray3(o, varStringl, fl, 
                                  nil, njl, nkl, cnl, eltTypel);
    if (resl == 1)
    {
        structVarString.push_back(varStringl);
        structF.push_back(fl);
        nit.push_back(nil); njt.push_back(njl), nkt.push_back(nkl);
        objs.push_back(o);
    }
    else if (resl == 2)
    {
        unstrVarString.push_back(varStringl);
        unstrF.push_back(fl);
        eltType.push_back(eltTypel); cnt.push_back(cnl);
        obju.push_back(o);
    }
    else printf("Warning: display: array " SF_D_ " is invalid. Discarded.\n", i);
  }

  Data* d = Data::getInstance();
