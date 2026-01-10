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
# include <stdexcept>
# include <cstring>
# include "GL/glew.h"
# include "ShaderUtil.h"
# include "ShaderObject.h"
using namespace CPlot;

//=============================================================================
ShaderObject::ShaderObject() :
  _shaderId(0),
  _nbReference(new unsigned int(1))
{}

//==============================================================================
ShaderObject::ShaderObject(const ShaderObject& obj) :
  _shaderId(obj._shaderId),
  _nbReference(obj._nbReference)
  
{
  *_nbReference += 1;
}
//==============================================================================
void ShaderObject::destroy()
{
  if (*_nbReference == 1)
  {
    delete _nbReference;
    glDeleteObjectARB(_shaderId);
    CHECK_GL_ERROR();
  }
  else
    *_nbReference -= 1;
}
//==============================================================================
ShaderObject::~ShaderObject()
{
  destroy();
}
//==============================================================================
const ShaderObject& ShaderObject::operator = (const ShaderObject& obj)
{
  if (this != &obj)
  {
    destroy();
    _nbReference = obj._nbReference;
    *_nbReference += 1;
    _shaderId = obj._shaderId;
  }
  return *this;
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
bool ShaderObject::compile(const char* src)
{
  bool is_compiled(false);
  GLint compiled = 0;
  GLint lgth = (GLint)strlen((const char*)src);
  glShaderSourceARB(_shaderId, 1, (const GLcharARB**)&src, &lgth);
  CHECK_GL_ERROR();

  glCompileShaderARB(_shaderId);
  CHECK_GL_ERROR();

  glGetObjectParameterivARB(_shaderId, GL_COMPILE_STATUS, &compiled);
  CHECK_GL_ERROR();
  
  if (compiled) is_compiled = true;

  return is_compiled;
}
//==============================================================================
void ShaderObject::load(const char* fileName)
{
  FILE* ptrFile = fopen(fileName, "r");
  if (ptrFile == NULL)
    throw std::runtime_error("File not found!");
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

  if (!compile(shadersource))
    throw std::runtime_error("Failed to compile shader source!");
  delete [] shadersource;
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

  glGetObjectParameterivARB(_shaderId, GL_OBJECT_COMPILE_STATUS_ARB, &isOK);
  CHECK_GL_ERROR();

  char* compiler_log = NULL;
  if (isOK != 1)
  {
    int infoLgth;
    glGetObjectParameterivARB(_shaderId, GL_OBJECT_INFO_LOG_LENGTH_ARB, &infoLgth);
    CHECK_GL_ERROR();
    if (infoLgth > 0) compiler_log = new char[infoLgth];
    if (compiler_log == NULL)
    {
      printf("ERROR: Could not allocate compiler_log buffer.\n");
      return errorMsg[2];
    }
    int msgLgth;
    glGetInfoLogARB(_shaderId, infoLgth, &msgLgth, compiler_log);
    CHECK_GL_ERROR();
    
  }
  std::string errMsg;
  if (compiler_log != NULL) errMsg = std::string(compiler_log);
  else errMsg = std::string(errorMsg[4]);
  delete [] compiler_log;
  return errMsg;
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
