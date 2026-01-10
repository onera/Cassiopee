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
#include <stdexcept>
#include "TesselationShaderManager.hpp"
#include "TesselationControlShader.hpp"
#include "TesselationEvaluationShader.hpp"
#include <string.h>
#include "../Data.h"
using namespace CPlot;

//=============================================================================
TesselationShaderManager::TesselationShaderManager() : 
  _shaderList(), 
  _currentActiveShader(0)
{}

//=============================================================================
TesselationShaderManager::~TesselationShaderManager()
{
  deactivate();
}

//=============================================================================
bool 
TesselationShaderManager::addFromFile( const char* tesselationBaseNameFile )
{
  std::string tesselationControlFileName   = std::string(tesselationBaseNameFile) + ".tcs";
  std::string tesselationEvaluationFileName= std::string(tesselationBaseNameFile) + ".tes";
  std::pair<std::shared_ptr<TesselationControlShader>,
            std::shared_ptr<TesselationEvaluationShader>> tesselation_shader;
  tesselation_shader.first.reset();
  tesselation_shader.second.reset();
  try 
  {
    tesselation_shader.first = 
       std::make_shared<TesselationControlShader>(ShaderObject::load(tesselationControlFileName.c_str()));
    tesselation_shader.second = 
       std::make_shared<TesselationEvaluationShader>(ShaderObject::load(tesselationEvaluationFileName.c_str()));
  }
  catch(std::runtime_error& err)
  {
    std::cerr << "Fail to add new tesselation shaders: " << tesselationBaseNameFile
              << std::endl;
    std::cerr << "\tCompilation error : " << err.what() << std::endl;
    if ( tesselation_shader.first != nullptr )
    {
      std::cerr << "\t Tesselation Control shader log:"
                << tesselation_shader.first->getCompilerLog() << std::endl;
    }
    if ( tesselation_shader.second != nullptr )
    {
      std::cerr << "\t Tesselation Evaluation shader log: "
                << tesselation_shader.second->getCompilerLog() << std::endl;
    }
    return false;
  }
  _shaderList.push_back(tesselation_shader);
  return true;
}
//=============================================================================
int TesselationShaderManager::load()
{
  char tes[256*8];
  Data* d = Data::getInstance();
  char* path = d->ptrState->shaderPath;

  // - 1 - Triangle HO polynomial IP_{2} - 6 vertices
  strcpy(tes,path); strcat(tes, "triangle_p2");
  addFromFile(tes);
  // - 2 - Quadrangle HO polynomial - 8 or 9 vertices
  strcpy(tes,path); strcat(tes, "quadrangle_ho");
  addFromFile(tes);

  // - 3 - HO Line polynomial
  strcpy(tes,path); strcat(tes, "mesh_line_ho");
  addFromFile(tes);
  return 1;
}
