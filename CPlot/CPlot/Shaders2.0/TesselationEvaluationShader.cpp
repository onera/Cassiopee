#include "GL/glew.h"
#include "TesselationEvaluationShader.hpp"
#include "ShaderUtil.h"
# include <iostream>
# include <stdexcept>
using namespace CPlot;

//==============================================================================
TesselationEvaluationShader::TesselationEvaluationShader(): ShaderObject()
{
  _shaderId = glCreateShader(GL_TESS_EVALUATION_SHADER);
  CHECK_GL_ERROR();
}
//==============================================================================
TesselationEvaluationShader::TesselationEvaluationShader(const std::string& src):
  ShaderObject()
{
  _shaderId = glCreateShader(GL_TESS_EVALUATION_SHADER);
  bool success = this->compile(src);
  if (not success)
  {
    std::string error("Failed to compile tesselation evaluation shader :\n");
    error += src; 
    throw std::runtime_error(error.c_str());  
  }
}
//==============================================================================
TesselationEvaluationShader::~TesselationEvaluationShader()
{}
