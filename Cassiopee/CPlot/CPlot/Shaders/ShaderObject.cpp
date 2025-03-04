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
# include <stdexcept>
# include <cstring>
# include "GL/glew.h"
# include "ShaderUtil.h"
# include "ShaderObject.h"
# include <iostream>
# include <cassert>
# include <vector>
using namespace CPlot;

//=============================================================================
ShaderObject::ShaderObject() :
  _shaderId(0)
{}
//==============================================================================
ShaderObject::~ShaderObject()
{
  glDeleteObjectARB(_shaderId);
  CHECK_GL_ERROR();
}
//==============================================================================
unsigned E_LONG ShaderObject::getFileLength(FILE* ptrFile)
{ 
  fseek(ptrFile, 0L, SEEK_END);
  unsigned E_LONG len = ftell(ptrFile);
  fseek(ptrFile, 0L, SEEK_SET);
  unsigned E_LONG pos = ftell(ptrFile);
  return (unsigned E_LONG)(len-pos);
}

//==============================================================================
bool ShaderObject::compile(const std::string& src)
{
  bool is_compiled(false);
  GLint compiled = 0;
  GLint lgth = (GLint)src.size();
  const GLchar* pt_src = (const GLchar*)src.data();
  glShaderSourceARB(_shaderId, 1, &pt_src, &lgth);
  CHECK_GL_ERROR();

  glCompileShaderARB(_shaderId);
  CHECK_GL_ERROR();

  glGetObjectParameterivARB(_shaderId, GL_COMPILE_STATUS, &compiled);
  CHECK_GL_ERROR();
  
  if (compiled == GL_TRUE) is_compiled = true;
  /*else
  {
    GLint maxLength = 0;
    glGetShaderiv(_shaderId, GL_INFO_LOG_LENGTH, &maxLength);

    // The maxLength includes the NULL character
    std::cout << "Error in shader :\n" << src << std::endl;
    std::vector<GLchar> errorLog(maxLength);
    glGetShaderInfoLog(_shaderId, maxLength, &maxLength, errorLog.data());
    for ( const auto& str : errorLog )
      std::cout << str;
    std::cout << std::endl;
  }*/

  return is_compiled;
}
//==============================================================================
std::string ShaderObject::load(const char* fileName)
{
  FILE* ptrFile = fopen(fileName, "r");
  if (ptrFile == NULL)
  {
    char errorMsg[1024];
    sprintf(errorMsg, "File %s not found !", fileName);
    throw std::runtime_error(std::string(errorMsg));
  }
  unsigned E_LONG lgth = getFileLength(ptrFile);
  if (lgth == 0) throw std::runtime_error("Empty Shader file!");
  char* shadersource = new char[lgth+1];
  if (shadersource == NULL) 
    throw std::runtime_error("Unable to reserve memory for shader source!");
  char c = fgetc(ptrFile);
  int i = 0;
  while (c != EOF) { shadersource[i] = c; i++; c = fgetc(ptrFile); }
  shadersource[i] = '\0';
  shadersource[lgth] = 0;// Max value for end of file...
  std::string s_source(shadersource);
  /*if (!compile(shadersource))
    throw std::runtime_error("Failed to compile shader source!");*/
  delete [] shadersource;
  return s_source;
}
//==============================================================================
static const char* errorMsg[] = {
  "Not a valid program object !",
  "Not a valid object !",
  "Out of memory",
  "Unknown compiler error !",
  "Compiler log is not available !",  
  NULL
};

//=============================================================================
std::string ShaderObject::getCompilerLog(void) const
{
  GLint isOK = 0;	

  if (_shaderId == 0) return std::string(errorMsg[1]);

  glGetShaderiv( _shaderId, GL_COMPILE_STATUS, &isOK );
  CHECK_GL_ERROR();

  if (isOK == GL_FALSE )
  {
    GLint infoLgth;
    glGetShaderiv( _shaderId, GL_INFO_LOG_LENGTH, &infoLgth );
    std::string slog(infoLgth+1,'\0' );
    glGetShaderInfoLog(_shaderId, infoLgth, &infoLgth, &slog[0]);
    CHECK_GL_ERROR();
    return slog;    
  }
  return std::string("No log available.");
}
//==============================================================================
unsigned int
ShaderObject::getAttributeLocation(const char* attributeName) const
{
  GLint idVar = glGetAttribLocationARB(_shaderId, attributeName);
  return (unsigned int)idVar;
}
//==============================================================================
unsigned int
ShaderObject::getUniformLocation(const char* uniformName) const
{
  GLint idVar = glGetUniformLocationARB(_shaderId, uniformName);
  return (unsigned int)idVar;
}
