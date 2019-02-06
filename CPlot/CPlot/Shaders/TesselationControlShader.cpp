#include "GL/glew.h"
#include "TesselationControlShader.hpp"
#include "ShaderUtil.h"
# include <iostream>
# include <stdexcept>
using namespace CPlot;

//==============================================================================
TesselationControlShader::TesselationControlShader(): ShaderObject()
{
  _shaderId = glCreateShader(GL_TESS_CONTROL_SHADER);
  CHECK_GL_ERROR();
}
//==============================================================================
TesselationControlShader::TesselationControlShader(const std::string& src ):
  ShaderObject()
{
  _shaderId = glCreateShader(GL_TESS_CONTROL_SHADER);
  bool success = this->compile(src);
  if (not success)
  {
    std::string error("Failed to compile tesselation control shader :\n");
    error += src; 
    throw std::runtime_error(error.c_str());  
  }
}
//==============================================================================
TesselationControlShader::~TesselationControlShader()
{}
